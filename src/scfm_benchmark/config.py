"""Single source of truth for environment paths and run parameters.

Every script or notebook that needs the project root, the data root, the logs
dir, or any tunable parameter should resolve it through this module instead of
hardcoding an absolute path. This is the one place the base path is defined.

Resolution order (lowest to highest precedence):
    1. configs/defaults.yaml            (environment-agnostic params)
    2. configs/<env>.yaml               (hpcc | local | cloud overlay)
    3. environment variables            (PROJECT_ROOT / SCFM_DATA_ROOT / SCFM_LOGS_DIR)

The active environment is chosen by the SCFM_ENV env var; if unset it is
auto-detected: `hpcc` when /bigdata exists (the UCR cluster), else `local`.

Usage
-----
    from scfm_benchmark.config import paths, settings
    df.to_csv(paths.data_root / "processed" / "x.csv")
    seed = settings["random_seed"]

Migration note: the numbered scripts under scripts/NN_*/ still carry their own
hardcoded /bigdata absolutes as a fallback (they run on HPCC where that path is
real). Migrate a script to `paths` opportunistically the next time it is edited;
see docs/decisions.md (2026-06-06 path-centralization decision).
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

# REPO_ROOT = the git checkout root (src/scfm_benchmark/config.py -> parents[2]).
REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIGS_DIR = REPO_ROOT / "configs"

_VALID_ENVS = ("hpcc", "local", "cloud")


def active_env() -> str:
    """Name of the active environment overlay (hpcc | local | cloud)."""
    env = os.environ.get("SCFM_ENV", "").strip().lower()
    if env in _VALID_ENVS:
        return env
    if env:
        raise ValueError(f"SCFM_ENV={env!r} not in {_VALID_ENVS}")
    # Auto-detect: the UCR cluster mounts /bigdata; a laptop does not.
    return "hpcc" if Path("/bigdata").exists() else "local"


def _load_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return yaml.safe_load(path.read_text()) or {}


def _deep_merge(base: dict[str, Any], over: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge `over` onto `base` (over wins on conflicts)."""
    out = dict(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def load_settings(env: str | None = None) -> dict[str, Any]:
    """Merged config dict: defaults.yaml overlaid by the env config."""
    env = env or active_env()
    settings = _deep_merge(
        _load_yaml(CONFIGS_DIR / "defaults.yaml"),
        _load_yaml(CONFIGS_DIR / f"{env}.yaml"),
    )
    settings["_env"] = env
    return settings


@dataclass(frozen=True)
class Paths:
    """Resolved filesystem roots for the active environment."""

    env: str
    project_root: Path
    data_root: Path
    logs_dir: Path

    def resolve(self, *parts: str) -> Path:
        """project_root / parts — the canonical way to build a project path."""
        return self.project_root.joinpath(*parts)


def load_paths(env: str | None = None) -> Paths:
    """Resolve project/data/logs roots, honoring env-var overrides."""
    s = load_settings(env)
    env = s["_env"]
    project_root = Path(os.environ.get("PROJECT_ROOT", s.get("project_root", str(REPO_ROOT))))
    data_root = Path(os.environ.get("SCFM_DATA_ROOT", s.get("data_root", str(project_root / "data"))))
    logs_dir = Path(os.environ.get("SCFM_LOGS_DIR", s.get("logs_dir", str(project_root / "logs"))))
    return Paths(env=env, project_root=project_root, data_root=data_root, logs_dir=logs_dir)


# Module-level singletons for the common case `from scfm_benchmark.config import paths`.
settings = load_settings()
paths = load_paths()


if __name__ == "__main__":  # quick introspection: `python -m scfm_benchmark.config`
    print(f"env         : {paths.env}")
    print(f"project_root: {paths.project_root}")
    print(f"data_root   : {paths.data_root}")
    print(f"logs_dir    : {paths.logs_dir}")
