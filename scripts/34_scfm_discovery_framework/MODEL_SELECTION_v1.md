# scFM backbone selection for the discovery framework — limited shortlist v1

Prepared per `review:20261012_140650` item 2-3: a bounded shortlist (≤2-3
candidates) with verified checkpoint/species/interpretation/engineering status, and
one recommended primary + one necessary alternative control. **Not a final decision —
the specific model/checkpoint is the owner's to confirm** (per that review's own
instruction: "把具体模型/checkpoint 推荐交给用户确认；不要把'单模型方向已批准'误当作某个模型已经选定").

No new SLURM job was run to produce this document. Every claim below is either (a) an
official source fetched directly this round (quoted, not relayed from the prior
review), or (b) a fact read directly from this repo's own files or from read-only
`ls`/`find`/log-tail commands on HPCC this round — never from the old README's
"production-ready" framing (per `review:20261012_133205` §5's caution about
that).

**Revised in place (2026-10-12) per `review:20261012_145004` §3.** That review
found the first version conflated scPRINT (v1) and scPRINT-2 evidence, understated
UCE's own species evidence, excluded the whole TranscriptFormer family on the basis
of one historical human-only run, and overstated how close the logged mouse attempts
came to a working result. All four are corrected below; the "most advanced attempt"/
"closest to succeeding" framing is retracted, not merely softened — the review is
correct that no evidence establishes either candidate as nearer to working than the
other, and a progress bar, an accuracy metric, or a "saving" print line is not proof
that the embedding was computed correctly end-to-end.

**Revised again (2026-10-12) per `review:20261012_153316`.** That review
closed scPRINT-2's "open item" by extracting the checkpoint FILE's own pickled
hyperparameters (not prose about it) and finding an `organisms` list there.
Independently re-verified this round with a separate method (a raw byte scan of
`archive/data.pkl` for `NCBITaxon:*` patterns — no pickle deserialization, no
REDUCE/GLOBAL/BUILD opcode executed): exactly 16 unique NCBITaxon IDs are present,
including `NCBITaxon:9606` (human) and `NCBITaxon:10090` (mouse), matching the
pickle disassembly's own 16-element `organisms` list exactly. scPRINT-2's
human/mouse species status is now **`checkpoint_metadata_confirmed`** (a new,
stronger evidence tier than README-level `documented_confirmed` — see
`model_registry.py`) — this is checkpoint metadata evidence that both species are
in the declared coverage list, not a per-species training-cell count, train/val
split, or any downstream efficacy claim, none of which this changes. Also
corrected: the required input gene-ID namespace is **Ensembl ID or HUGO symbol**
(re-verified this round by re-fetching the README FAQ directly), not `gene_symbol`
only as the first version wrongly required. That review also reports reading the
UCE Nature paper (`s41586-026-10689-z`) directly and finding its species paragraph
names human and mouse explicitly; my own two fetch attempts this round were both
blocked (an auth redirect to `idp.nature.com`) and a bioRxiv retry was rate-limited
(HTTP 429) — I could not independently reproduce that specific read, so UCE's
species evidence below is left at `documented_unconfirmed` from my own verification
standard, with this attribution noted rather than silently adopted.

## Verification status

UNVERIFIED (recommendation only; no checkpoint pinned, no adapter built, no new
inference run). Status label follows the three states the review asked for —
**documented** (official source says so), **locally available** (the actual weights/
env exist on this HPCC account right now), **tested** (an actual run against this
project's own data produced a verifiable, readable artifact) — plus a fourth state
this investigation surfaced and that matters for the mouse decision: **attempted,
not yet working** (a real run was tried against real data and did not produce a
verifiable artifact). Absence of a success is not evidence of impossibility (DEVELOPMENT.md
rule 3) — each "attempted, not yet working" row below states the actual failure
point, not a verdict on the model.

## Why species is the actual decision axis here

This project's *old* 5-model benchmark (STATE/UCE/Geneformer/scPRINT-2/scGPT, plus
STACK/TranscriptFormer-Sapiens) already generated real embeddings against the real
human platelet/COVID canonical data for every one of those models —
`/bigdata/.../02_EMBEDDINGS/option_c/{state,uce,geneformer,scprint2,scgpt,stack,
transcriptformer}_covid_canonical.h5ad` all exist on disk, all multi-GB, dated
2026-05-01 through 06-08. That these files exist and are the right size for real
data is a provenance fact, not a validation of them against the *new* framework's
monocyte input, cell-ID alignment, or current QC version — none of that has been
checked this round. So for the narrower purpose of this document (choosing which
backbone is *structurally* able to serve mouse at all), the human side is not the
differentiator, because something has already run against human data in this repo
for essentially every candidate. The new, genuinely open question is **mouse**
support, because the framework's stated order (`review:20261012_133528`) is
one primary backbone with documented species coverage. Model selection therefore turns on cross-species evidence,
not on human performance, and not on which model's old human artifacts happen to be
on disk.

## Candidates considered and ruled out structurally (human-only)

| Model | Official source (fetched this round) | Verdict |
|---|---|---|
| **Geneformer** | `huggingface.co/ctheodoris/Geneformer`: "pretrained ... on a large-scale corpus of **human** single cell transcriptomes" (V1: ~30M human cells; V2: ~104M human, non-cancer). No mention anywhere on the card of a mouse-trained checkpoint or mouse vocabulary. | **Human-only by the official checkpoint.** Already has real, tested human embeddings in this repo (`geneformer_covid_canonical.h5ad`, 600MB) and a dedicated `geneformer_py310` conda env — cheap to keep as an internal human reference point, but structurally disqualified as the one backbone that also has to serve mouse. |
| **scGPT** | `github.com/bowang-lab/scGPT` README, Model Zoo: only the **whole-human** checkpoint ("Pretrained on 33 million normal human cells") names a species explicitly; the tissue-specific checkpoints (brain/blood/heart/lung/kidney/pan-cancer) don't state species but give no indication of mouse. No checkpoint names mouse anywhere on the page. | **No official mouse checkpoint found.** Real tested human embeddings exist (`scgpt_covid_canonical.h5ad`, 2.6GB) and two dedicated envs (`scgpt`, `scgpt_finetune`) — same disposition as Geneformer: good human reference, disqualified as the cross-species backbone. |
| **STATE** (Arc Institute) | Web search of Arc's own announcement pages: training described as "167 million cells" observational + "100 million" perturbational data, explicitly across "**70 human** cell contexts." The only mouse mention found is a third-party LinkedIn post describing *downstream application* to mouse experiments, not training data. | **No official mouse-training claim found.** Has a dedicated `state_env` and the largest existing human embedding file in this repo (`state_covid_canonical.h5ad`, 5.4GB) — strong human reference, but there is no official basis to claim mouse support, so it is excluded from the cross-species decision on the same grounds as the two above. |

These three are not weaker models — they are simply not candidates for the *one*
backbone that must also cover mouse. Their existing human embeddings remain
available as internal reference points if the paper wants one — using them would
still require re-verifying them against the new monocyte input/QC version, not a
zero-cost drop-in, since none of that re-verification has been done.

## The genuine cross-species candidates

### scPRINT-2 — recommended primary

**scPRINT and scPRINT-2 are different generations and must not be cited
interchangeably.** The first version of this document quoted the *scPRINT (v1)*
README's species sentence while recommending *scPRINT-2*, and listed `large-v1.ckpt`
(a scPRINT v1 weight) and `small-v2.ckpt` (a scPRINT-2 weight) as if they were
interchangeable options for one model. They are not. Everything below is re-sourced
from the dedicated `cantinilab/scPRINT-2` repository and the local `small-v2.ckpt`
checkpoint specifically.

**Documented (README):** the scPRINT-2 README (fetched directly this round) states
it is *"pretrained on more than 350 million cells from more than 22,000 datasets
and 16 species"* — but the README itself does not name which 16 species; it tells
the reader to "check in the model.organisms or in our manuscript." A guessed repo
path for `model.organisms` 404'd, and the manuscript (bioRxiv) was rate-limited
twice. **This was correctly logged as an open item in the prior revision of this
document.**

**Resolved this round, via the checkpoint file itself, not the README.**
`review:20261012_153316` located the actual `organisms` list inside
`small-v2.ckpt`'s own saved hyperparameters (read by the official
`cantinilab/scPRINT-2` model source, which loads this field from the checkpoint at
load time — so it is a real model/checkpoint attribute, not a guessable separate
file). Independently re-verified this round with a different, equally non-executing
method — a raw byte scan of the checkpoint's `archive/data.pkl` for `NCBITaxon:*`
substrings, with no pickle deserialization and no REDUCE/GLOBAL/BUILD opcode ever
run: **exactly 16 unique NCBITaxon IDs are present, including `NCBITaxon:9606`
(human, 1 occurrence) and `NCBITaxon:10090` (mouse, 1 occurrence)**, matching the
pickle disassembly's own 16-element `organisms` list length exactly. **scPRINT-2's
human/mouse species status is now `checkpoint_metadata_confirmed`** — stronger than
README-level `documented_confirmed` because it is the specific checkpoint's own
declared coverage, not prose about the model family. This still does NOT establish
per-species training-cell counts, train/validation split, downstream
training representation, or any downstream efficacy — none of those were checked
and none follow from a species appearing in a coverage list.

What *is* solidly documented for scPRINT-2: species is selected via a checkable
ontology ID at inference (`--species "NCBITaxon:9606"` in the `embed`/`gninfer`
examples), not a free string; and the model's own stated purpose is gene-level, not
only cell-level — *"generate a gene network from any cell or cell cluster"* via a
dedicated `scprint2 gninfer` CLI verb, with a companion `GRnnData` package for
working with the resulting networks. This interface shape is the real reason
scPRINT-2 fits the framework's "cell state → gene program" step better than a
cell-embedding-only model — **a network/attention-derived gene connection is not
itself a causal regulatory relationship**; it is model-inferred structure that would
still need independent support before being read as biology.

**Checkpoint identity, verified this round:** the local cached file
`/bigdata/godziklab/shared/Xinru/.huggingface/scprint2_ckpts/small-v2.ckpt`
(889,706,878 bytes) hashes to
`sha256:2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822` —
**this matches exactly** the SHA-256 published on the official weight page
(`huggingface.co/jkobject/scPRINT/blob/main/small-v2.ckpt`, fetched this round,
same hash, "890 MB"). This is a real, checked identity binding, not an assumption.

**Locally available:** dedicated conda env `scprint2_env`, with `small-v2.ckpt`
cached as above.

**Tested (human):** real embeddings exist for this project's own COVID canonical
data (`scprint2_covid_canonical.h5ad`, 400MB) and were carried through the existing
(pre-framework) GSEA step — `05_RESULTS/option_c_gsea_fullrank/
gsea_scprint2_dim*_rank*_{KEGG_2021_Human,Reactome_2022,GO_Biological_Process_2023}.csv`
exist for 5 scPRINT-2 dimensions. This is a 2026-05/06 result against the old
platelet/COVID cohort and QC version, on this exact project's data — it is **not**
evidence that the new monocyte input, cell-ID alignment, or current QC version
reproduces anything; it establishes only that the scPRINT-2 environment and
inference path have worked in this repo before, nothing about the new framework's
specific inputs.

### UCE — recommended alternative/control

**Documented:** `github.com/snap-stanford/UCE` README: `--species` is a free-form
string (no `choices=` restriction) — species coverage lives in per-species asset
files, not a hardcoded list, so the *interface* does not itself name species. The
primary training-data evidence is the UCE paper itself (Nature, `s41586-026-10689-z`
— I hit an authentication redirect fetching it directly this round, so could not
quote it verbatim) and its earlier bioRxiv preprint (`10.1101/2023.11.28.568918`,
rate-limited (HTTP 429) on two direct fetch attempts this round). The best
corroboration I could independently obtain this round is a secondary source
describing the preprint's own dataset breakdown: **human and mouse cells together
account for 33.9 million of the paper's 36.2 million total training cells** — a far
more specific and substantial claim than "listed among 8 species," and not merely an
"open-ended list," correcting the first version of this document's framing. I was
not able to fetch the primary-source sentence verbatim this round, so this remains
secondary-sourced, not a direct quote.

There are two incompatible checkpoint architectures (4-layer, default; 33-layer,
separate download) whose embeddings the README explicitly warns "are not compatible
with" each other — the model used in this project's local install and in all mouse
attempts below is the 4-layer checkpoint.

**TranscriptFormer is a related, documented cross-species family that the first
version of this document did not carry into the final table** (it appeared only in
the February-round research, not here). Per `github.com/czi-ai/transcriptformer`
README (re-confirmed this round): **TF-Sapiens is human-only**; **TF-Exemplar**
("human plus four model organisms: mouse (M. musculus), zebrafish, fruit fly, and
C. elegans") and **TF-Metazoa** ("112 million cells spanning all twelve species")
both explicitly include mouse. This project has only ever run **TF-Sapiens**
(human-only) successfully (`covid_transcriptformer.h5ad`, per
`SESSION_SUMMARY.md`); a prior TF-Metazoa download attempt failed outright (0-byte
files, wrong CLI casing: `'TF-Metazoa'` vs. the required lowercase `'tf-metazoa'`).
**That TF-Sapiens-only history is not grounds to exclude the TF-Exemplar/TF-Metazoa
checkpoints from future consideration** — they are a real, documented, not-yet-tried
cross-species option, noted here rather than built out, per the review's own
"补一行说明即可，无须扩大为多模型实现或全面比较" instruction.

**Locally available:** this project's own UCE install
(`/bigdata/.../03_MODELS/uce/UCE/model_files/`) already has the mouse-specific
protein-embedding file (`Mus_musculus.GRCm39.gene_symbol_to_embedding_ESM2.pt`,
confirmed present by `ls` this round) — a mouse run does not need any new download.

**Tested (human):** real embeddings exist (`uce_covid_canonical.h5ad`, 5.0GB) and
downstream GSEA results exist for UCE dimensions too
(`gsea_uce_dim*_rank*_{KEGG,Reactome,GO}.csv`) — same caveat as scPRINT-2 above: a
2026 result on the old cohort/QC version, not evidence for the new framework's
inputs.

## Recommendation

**Primary: scPRINT-2, as a candidate worth a minimal validation first — not as an
already-proven choice.** scPRINT-2's gene-network-inference interface (`gninfer`,
`GRnnData`) is the closest documented match to the framework's own "cell state →
gene program" step of anything checked; its exact `small-v2.ckpt` checkpoint
identity is verified (local SHA-256 matches the official page exactly); and its
own checkpoint metadata now confirms human+mouse are in its declared 16-species
coverage (resolved this round, see above — `checkpoint_metadata_confirmed`).
**This is still not a claim that scPRINT-2 is proven superior to UCE, that its
downstream downstream efficacy is established, or that its logged mouse
attempt succeeded** — none of those are established, and a species appearing in a
checkpoint's coverage list says nothing about per-species training depth or actual
embedding quality on this project's specific mouse data.

**Alternative/control: UCE.** Grounds: the review's own ask for "一个必要时可用的代表对照"
(one available control if needed) — UCE has real tested human embeddings, real
downstream GSEA output, and the needed mouse protein-embedding asset already on
disk. Its species evidence remains at `documented_unconfirmed` from my own
verification standard (the review reports reading the Nature paper directly and
finding human+mouse named explicitly; I could not reproduce that specific read
myself this round — see above) — offered as a parallel option regardless, not a
weaker one.

**This is a recommendation, not a selection.** Per the review's own instruction,
pinning the exact model AND checkpoint, and deciding whether to debug an existing
failed mouse attempt or start fresh, is the owner's decision once the open species
question above is resolved — not something this document or any AI review decides.

## What this document does not establish

- Does not run any new inference, debug the scPRINT-2/UCE mouse attempts, or build
  any adapter. All evidence above is read-only (official-source fetches + `ls`/`find`/
  log-tail on existing files and already-completed jobs).
- scPRINT-2's checkpoint declares human+mouse coverage (`checkpoint_metadata_confirmed`,
  resolved this round); UCE's does not have an equivalent checkpoint-level check
  performed and remains `documented_unconfirmed` from my own verification. Neither
  establishes per-species training depth, train/val composition, or downstream
  downstream efficacy — this document does not resolve those.
- Does not establish that either scPRINT-2's or UCE's logged mouse attempt came
  closer to a working result than the other's — both are "historical attempt,
  result unverified," corrected down from the first version's "closest to
  succeeding" framing, which this review correctly identified as unsupported.
- Does not evaluate scPRINT-2/UCE against any biological ground truth, and does not
  treat a gene-network/attention-derived connection as a causal regulatory
  relationship — "fits the interpretation goal" and "best-documented species
  contract" are engineering/provenance judgments, not claims about biological
  effect (per the review's own caution: "'容易部署' '适合这个任务' 和'生物学效果更好'是不同结论").
- Does not treat the existing human-side canonical embeddings/GSEA output as
  validated for the new framework's monocyte input, cell-ID alignment, or current
  QC version — they are historical artifacts from the old cohort/QC version, real
  but not yet re-verified against anything the new framework will actually run.
- Does not address STATE/scPRINT-2/UCE/scGPT's human-side GSEA results' use of a
  live Enrichr API by library name (not a cached/versioned gene-set file) — this is a
  real gap against the new framework's own requirement
  (`review:20261012_133205` §3) and belongs to the enrichment-provider
  implementation task, not to this selection document.
- Does not re-examine `GSM5102900_HC1`/`GSM5102901_HC2` Gate 2 work, which is
  unaffected by this scope change.

## Reusable components inventory (informs, does not replace, the foundation build)

Read-only inspection this round, for the next implementation increment:

- `src/scfm_benchmark/config.py`: an existing `paths`/`settings` resolver
  (`SCFM_ENV`-aware, `/bigdata` auto-detect) — directly reusable for the new
  framework's path contracts, not a rewrite target.
- `scripts/24_option_c_canonical_input/code/gsea_per_fm_dim.py`: a working
  `gseapy.prerank`-based signed-ranking GSEA runner — reusable as the *ranking/
  prerank* mechanic, but it calls Enrichr's live API by library name
  (`KEGG_2021_Human`, etc.) with no cached gene-set file, no content hash, and no
  offline replay — **needs fixing**, not directly reusable, against the new
  framework's cached/versioned-provider requirement.
- `src/scfm_benchmark/api.py`/`cli.py` (`embed` command): hardcoded to the old HPCC
  sbatch submitters for the fixed 5-model set; `--h5ad/--out` are accepted but not
  actually forwarded to the adapter call. **Old scientific objective, not directly
  reusable** for the new discovery framework's adapter interface, though the
  argument-parsing shell is a reasonable reference.
- `02_EMBEDDINGS/option_c/*_covid_canonical.h5ad`: real, already-generated human
  embeddings for 7 models (state/uce/geneformer/scprint2/scgpt/stack/
  transcriptformer), against the **old** platelet/COVID cohort and QC version, not
  the new monocyte input. Historical assets with their own provenance (job IDs,
  dates) — not validated for, and not to be re-claimed as, new-framework output
  without re-running against the new input/QC version.

## Next implementation step

Per `review:20261012_140650` item 4, the next increment (once the owner
confirms or redirects the primary/control choice above) is to implement exactly one
production adapter for the chosen model and wire the minimal state→gene→pathway
evidence chain, before any pathway-provider/Docker expansion. Until that choice is
confirmed, continuing the already-authorized, model-agnostic foundation work (dataset/
checkpoint-species contracts, explicit `validate`/`plan` entry points, artifact
provenance) does not need to wait, since it does not depend on which backbone is
finally chosen.
