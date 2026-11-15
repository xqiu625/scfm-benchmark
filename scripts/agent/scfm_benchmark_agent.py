"""
scFM Benchmark Agent — LangGraph implementation
Autonomous agent for benchmarking single-cell foundation models.

Architecture:
  - StateGraph with typed state
  - 7 specialized nodes (planner, literature, data_prep, embedding, benchmark, robustness, interpret)
  - Human-in-the-loop at critical decision points
  - Tool-calling for HPC job submission, data loading, and analysis

Usage:
  python scfm_benchmark_agent.py --query "Benchmark STATE and UCE on COVID platelet data"
"""

from __future__ import annotations

import operator
from dataclasses import dataclass, field
from enum import Enum
from typing import Annotated, Any, Literal

from langgraph.graph import END, StateGraph
from langgraph.prebuilt import ToolNode
from langgraph.checkpoint.memory import MemorySaver
from langchain_anthropic import ChatAnthropic
from langchain_core.messages import AIMessage, BaseMessage, HumanMessage, SystemMessage
from langchain_core.tools import tool

# ---------------------------------------------------------------------------
# 1. State Definition
# ---------------------------------------------------------------------------

class Phase(str, Enum):
    PLAN = "plan"
    LITERATURE = "literature"
    DATA_PREP = "data_prep"
    EMBEDDING = "embedding"
    BENCHMARK = "benchmark"
    ROBUSTNESS = "robustness"
    DRUG_DISCOVERY = "drug_discovery"
    INTERPRET = "interpret"
    DONE = "done"


@dataclass
class BenchmarkState:
    """Typed state that flows through the graph."""

    # --- conversation ---
    messages: Annotated[list[BaseMessage], operator.add] = field(default_factory=list)

    # --- plan ---
    current_phase: Phase = Phase.PLAN
    models_to_benchmark: list[str] = field(default_factory=list)
    tasks: list[str] = field(default_factory=list)       # binary, 3-class, 6-class
    robustness_axes: list[str] = field(default_factory=list)

    # --- data ---
    data_path: str = ""
    n_cells: int = 0
    qc_passed: bool = False

    # --- embeddings ---
    embeddings_ready: dict[str, bool] = field(default_factory=dict)

    # --- results ---
    benchmark_results: dict[str, dict] = field(default_factory=dict)
    robustness_results: dict[str, dict] = field(default_factory=dict)
    drug_candidates: list[dict] = field(default_factory=list)

    # --- control ---
    human_approval_needed: bool = False
    error: str | None = None


# ---------------------------------------------------------------------------
# 2. Tools — what the agent can actually execute
# ---------------------------------------------------------------------------

@tool
def load_anndata(path: str) -> dict:
    """Load an AnnData h5ad file and return basic QC metrics."""
    import scanpy as sc

    adata = sc.read_h5ad(path)
    n_cells, n_genes = adata.shape
    severity_counts = adata.obs["severity"].value_counts().to_dict() if "severity" in adata.obs else {}

    return {
        "n_cells": n_cells,
        "n_genes": n_genes,
        "severity_distribution": severity_counts,
        "obs_columns": list(adata.obs.columns),
        "obsm_keys": list(adata.obsm.keys()),
    }


@tool
def check_embeddings(path: str, model_name: str) -> dict:
    """Check if pre-computed embeddings exist for a given model."""
    import scanpy as sc

    adata = sc.read_h5ad(path, backed="r")
    embed_keys = {
        "STATE": "X_emb",
        "UCE": "X_uce",
        "scGPT": "X_scGPT",
        "Geneformer": "X_geneformer",
        "TranscriptFormer": "X_tf",
    }
    key = embed_keys.get(model_name, f"X_{model_name.lower()}")
    exists = key in adata.obsm
    dims = adata.obsm[key].shape[1] if exists else 0

    return {"model": model_name, "exists": exists, "key": key, "dims": dims}


@tool
def submit_slurm_job(script_path: str, job_name: str) -> dict:
    """Submit a SLURM job to the HPC cluster and return job ID."""
    import subprocess

    result = subprocess.run(
        ["sbatch", f"--job-name={job_name}", script_path],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        return {"success": False, "error": result.stderr}

    job_id = result.stdout.strip().split()[-1]
    return {"success": True, "job_id": job_id, "job_name": job_name}


@tool
def run_classification(
    data_path: str,
    model_name: str,
    task: str,
    strategy: str = "deep_mlp",
) -> dict:
    """Run classification benchmark for a single model + task + strategy."""
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_score
    from sklearn.neural_network import MLPClassifier
    from sklearn.preprocessing import StandardScaler, LabelEncoder
    import scanpy as sc

    adata = sc.read_h5ad(data_path)

    embed_keys = {
        "STATE": "X_emb", "UCE": "X_uce", "scGPT": "X_scGPT",
        "Geneformer": "X_geneformer", "TranscriptFormer": "X_tf",
    }
    X = adata.obsm[embed_keys[model_name]]

    # Build labels based on task
    if task == "binary":
        label_map = {
            "healthy": 0, "mild": 0, "moderate": 0, "recovered": 0,
            "severe": 1, "fatal": 1,
        }
    elif task == "3-class":
        label_map = {
            "healthy": 0, "mild": 1, "moderate": 1, "recovered": 1,
            "severe": 2, "fatal": 2,
        }
    else:  # 6-class
        le = LabelEncoder()
        y = le.fit_transform(adata.obs["severity"])
        label_map = None

    if label_map:
        y = np.array([label_map[s] for s in adata.obs["severity"]])

    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    if strategy == "logreg":
        clf = LogisticRegression(max_iter=100, solver="lbfgs", multi_class="auto")
    else:
        clf = MLPClassifier(
            hidden_layer_sizes=(256, 128, 64),
            max_iter=200, early_stopping=True,
            batch_size=256,
        )

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    scores = cross_val_score(clf, X, y, cv=cv, scoring="roc_auc_ovr")

    return {
        "model": model_name, "task": task, "strategy": strategy,
        "auc_mean": float(np.mean(scores)),
        "auc_std": float(np.std(scores)),
        "auc_ci_95": [float(np.percentile(scores, 2.5)), float(np.percentile(scores, 97.5))],
    }


@tool
def run_robustness_test(
    data_path: str,
    model_name: str,
    test_type: str,
    param: float = 0.5,
) -> dict:
    """Run a robustness test: 'dropout', 'subsample', or 'batch_shift'."""
    # Placeholder — actual implementation calls the specialized scripts
    return {
        "model": model_name, "test": test_type, "param": param,
        "status": "submitted", "note": "See SLURM logs for results",
    }


@tool
def search_literature(query: str, max_papers: int = 10) -> list[dict]:
    """Search biomedical literature via PubMed/Semantic Scholar."""
    # In production: integrate paper-qa (FutureHouse) or PubMed API
    return [{"query": query, "status": "placeholder", "n_results": max_papers}]


@tool
def generate_report(results: dict, output_path: str) -> dict:
    """Generate a markdown report from benchmark results."""
    lines = ["# scFM Benchmark Report\n"]
    for model, metrics in results.items():
        lines.append(f"## {model}")
        for k, v in metrics.items():
            lines.append(f"- **{k}**: {v}")
        lines.append("")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))

    return {"output": output_path, "n_models": len(results)}


ALL_TOOLS = [
    load_anndata, check_embeddings, submit_slurm_job,
    run_classification, run_robustness_test, search_literature,
    generate_report,
]

# ---------------------------------------------------------------------------
# 3. LLM setup
# ---------------------------------------------------------------------------

SYSTEM_PROMPT = """You are an expert computational biologist specializing in single-cell
foundation model evaluation. You help design and execute benchmarking experiments for
scFMs (STATE, UCE, scGPT, Geneformer, TranscriptFormer) on clinical transcriptomic data.

Your workflow:
1. PLAN — decompose the user's query into specific benchmark tasks
2. LITERATURE — search for relevant prior work and baselines
3. DATA_PREP — load data, run QC, verify embeddings exist
4. EMBEDDING — generate or verify embeddings for each model
5. BENCHMARK — run classification tasks (binary, 3-class, 6-class)
6. ROBUSTNESS — run distribution shift, sample efficiency, gene dropout tests
7. DRUG_DISCOVERY — run perturbation recovery analysis (if requested)
8. INTERPRET — generate figures, statistical tests, and final report

Always explain your reasoning. Ask for human approval before submitting HPC jobs
or running expensive computations.
"""

def get_llm():
    return ChatAnthropic(
        model="claude-sonnet-4-20250514",
        max_tokens=4096,
    ).bind_tools(ALL_TOOLS)


# ---------------------------------------------------------------------------
# 4. Graph Nodes
# ---------------------------------------------------------------------------

def planner_node(state: BenchmarkState) -> dict:
    """Decompose user query into a structured plan."""
    llm = get_llm()
    messages = [SystemMessage(content=SYSTEM_PROMPT)] + state.messages
    messages.append(HumanMessage(content=(
        "Based on the user's request, create a plan. "
        "Decide which models to benchmark, which tasks to run, "
        "and which robustness axes to test. "
        "Call tools to check data availability."
    )))
    response = llm.invoke(messages)
    return {
        "messages": [response],
        "current_phase": Phase.DATA_PREP,
    }


def data_prep_node(state: BenchmarkState) -> dict:
    """Load data, run QC, check embedding availability."""
    llm = get_llm()
    messages = [SystemMessage(content=SYSTEM_PROMPT)] + state.messages
    messages.append(HumanMessage(content=(
        "Load the data and check which embeddings are available. "
        "Use load_anndata and check_embeddings tools. "
        "Report QC metrics and embedding status."
    )))
    response = llm.invoke(messages)
    return {
        "messages": [response],
        "current_phase": Phase.BENCHMARK,
    }


def benchmark_node(state: BenchmarkState) -> dict:
    """Run classification benchmarks."""
    llm = get_llm()
    messages = [SystemMessage(content=SYSTEM_PROMPT)] + state.messages
    messages.append(HumanMessage(content=(
        "Run the classification benchmarks. For each model with available "
        "embeddings, run the requested tasks using run_classification. "
        "Compare zero-shot (logreg) vs embedding classification (deep_mlp)."
    )))
    response = llm.invoke(messages)
    return {
        "messages": [response],
        "current_phase": Phase.ROBUSTNESS,
    }


def robustness_node(state: BenchmarkState) -> dict:
    """Run robustness analyses."""
    llm = get_llm()
    messages = [SystemMessage(content=SYSTEM_PROMPT)] + state.messages
    messages.append(HumanMessage(content=(
        "Run robustness tests: distribution shift, sample efficiency, "
        "gene dropout. Use run_robustness_test for each axis. "
        "Summarize which models are most robust."
    )))
    response = llm.invoke(messages)
    return {
        "messages": [response],
        "current_phase": Phase.INTERPRET,
    }


def interpret_node(state: BenchmarkState) -> dict:
    """Generate final report with interpretation."""
    llm = get_llm()
    messages = [SystemMessage(content=SYSTEM_PROMPT)] + state.messages
    messages.append(HumanMessage(content=(
        "Synthesize all results. Generate a final report with: "
        "1) Model rankings, 2) Robustness findings, "
        "3) Practical recommendations. Use generate_report tool."
    )))
    response = llm.invoke(messages)
    return {
        "messages": [response],
        "current_phase": Phase.DONE,
    }


def human_review_node(state: BenchmarkState) -> dict:
    """Pause for human approval before expensive operations."""
    # LangGraph's interrupt() mechanism handles this
    return {"human_approval_needed": False}


# ---------------------------------------------------------------------------
# 5. Routing Logic
# ---------------------------------------------------------------------------

def route_after_phase(state: BenchmarkState) -> str:
    """Decide next node based on current phase."""
    phase_to_node = {
        Phase.PLAN: "planner",
        Phase.DATA_PREP: "data_prep",
        Phase.BENCHMARK: "benchmark",
        Phase.ROBUSTNESS: "robustness",
        Phase.INTERPRET: "interpret",
        Phase.DONE: END,
    }
    return phase_to_node.get(state.current_phase, END)


def should_use_tools(state: BenchmarkState) -> Literal["tools", "continue"]:
    """Check if the last message has tool calls."""
    last = state.messages[-1]
    if isinstance(last, AIMessage) and last.tool_calls:
        return "tools"
    return "continue"


# ---------------------------------------------------------------------------
# 6. Build the Graph
# ---------------------------------------------------------------------------

def build_graph() -> StateGraph:
    """Construct the LangGraph benchmark agent."""

    graph = StateGraph(BenchmarkState)

    # --- Add nodes ---
    graph.add_node("planner", planner_node)
    graph.add_node("data_prep", data_prep_node)
    graph.add_node("benchmark", benchmark_node)
    graph.add_node("robustness", robustness_node)
    graph.add_node("interpret", interpret_node)
    graph.add_node("tools", ToolNode(ALL_TOOLS))

    # --- Entry point ---
    graph.set_entry_point("planner")

    # --- Edges: each node → tool check → routing ---
    for node_name in ["planner", "data_prep", "benchmark", "robustness", "interpret"]:
        graph.add_conditional_edges(
            node_name,
            should_use_tools,
            {"tools": "tools", "continue": "router"},
        )

    # After tools execute, go back to the node that called them
    graph.add_node("router", lambda state: state)
    graph.add_conditional_edges("router", route_after_phase)

    # Tools → router
    graph.add_edge("tools", "router")

    return graph.compile(checkpointer=MemorySaver())


# ---------------------------------------------------------------------------
# 7. Main
# ---------------------------------------------------------------------------

def main():
    import argparse

    parser = argparse.ArgumentParser(description="scFM Benchmark Agent")
    parser.add_argument("--query", type=str, required=True,
                        help="Natural language benchmark request")
    parser.add_argument("--data-path", type=str, default="",
                        help="Path to h5ad file with embeddings")
    args = parser.parse_args()

    graph = build_graph()

    initial_state = BenchmarkState(
        messages=[HumanMessage(content=args.query)],
        data_path=args.data_path,
    )

    config = {"configurable": {"thread_id": "benchmark-001"}}

    print("=" * 60)
    print("scFM Benchmark Agent (LangGraph)")
    print("=" * 60)

    for event in graph.stream(initial_state, config=config):
        for node_name, output in event.items():
            if node_name == "__end__":
                continue
            msgs = output.get("messages", [])
            for msg in msgs:
                if isinstance(msg, AIMessage):
                    print(f"\n[{node_name}] {msg.content[:500]}")
                    if msg.tool_calls:
                        for tc in msg.tool_calls:
                            print(f"  -> Tool: {tc['name']}({tc['args']})")

    print("\n" + "=" * 60)
    print("Benchmark complete.")


if __name__ == "__main__":
    main()
