# LROT

This directory is an independent repository root. Frozen GO/MGI inputs are bundled
in `data/go/`, and Supplementary File S1 is bundled in `supplementary_files/`.
No sibling delivery directories or separate resource archive are required.
See [repository preparation and checks](docs/repository_upload.md).

Ligand–receptor-guided optimal transport for spatial transcriptomics section alignment.

LROT combines gene-expression similarity, a manually specified ligand–receptor
(LR) compatibility prior and within-section spatial structure to obtain a soft
transport plan between two sections. The package contains the implementation,
fixed execution protocol, selected saved results, figure inputs and verification
records. The LROT solver itself does not train a neural network.

[中文说明](README_zh.md) · [Running guide](RUN_REPRO.md) ·
[Figure commands](docs/figure_entrypoints.md) · [Verification scope](docs/release_readiness.md)

## Package structure

| Location, relative to the code directory | Contents |
|---|---|
| `paper_protocol/` | Source verification, explicit numerical binding and the locked execution wrapper |
| `provenance/numerical_validation_snapshot/code/` | Preserved numerical implementation used by the bound protocol |
| `experiments/` | Case-specific preparation, solving and evaluation entrypoints |
| `figures/` | Approved images, saved drawing inputs, manifests and supported renderers |
| `lrot_output/` | Selected numerical summaries, saved arrays and seed-pair records |
| `data/`, `cancer_data/`, `dlpfc_data/` | Input identities, checksums, preparation instructions and resource notices |
| `docs/`, `tools/` | Environment profiles, guides and package/input verification tools |
| `data/go/` | Fixed ontology and annotation inputs with their provenance and notices |
| `supplementary_files/` | Supplementary File S1 workbook |

The original tissue matrices and CellChatDB binaries are obtained separately.
Manuscript DOCX/PDF files, virtual environments, credentials, upload utilities,
temporary logs and private backups are not included. Third-party source links in
the data and licence documentation identify required resources; they are not a
hosting address for this package.

## Complete method workflow

### 1. Identify and prepare the inputs

For sections A and B, prepare ordered spot-by-gene expression arrays, matching
spot coordinates and a common ordered gene list. Sample identities, barcodes,
coordinate units, transformations and preprocessing must follow the selected
case's input binding. Labels are required only for annotation-based evaluation.

Raw counts, transformed expression and retained plotting arrays are different
inputs. Use the relevant loader and recorded configuration; do not apply a new
normalization, reorder sampled indices or change precision to match a result.
Keep each recorded seed, sample direction and gene order fixed.

The [input guide](data/README.md) and [file manifest](data/input_files.json)
provide source pointers and checksums. DLPFC cases are identified by section IDs;
breast-cancer cases use Block A Sections 1 and 2. The mouse perturbation case
uses one measured expression dataset with generated coordinates and imposed
transformations, rather than two independently measured serial sections.
Simulation cases follow the fixed design and registered seed-pair file; here
"registered" means recorded, not preregistered.

### 2. Select the species-specific LR prior

Use the existing mouse or human pair table and its actual weights. These are
manually assigned compatibility priors, not measured interaction strengths or
database confidence scores. Missing genes are skipped as implemented; gene
symbol handling and feature selection remain case-specific. Do not substitute
another species table or infer orthology from capitalization alone.

Pair-level interpretation and resource limitations are documented in
[the LR resource guide](docs/lr_resources.md) and
[the evidence record](data/lr_resource_evidence.json). Supplementary File S1 is
the [mouse/human workbook](supplementary_files/LROT_Supplementary_Tables_S11_S12.xlsx).

### 3. Construct the expression, LR and spatial channels

The bound numerical source defines the following operations:

* **Expression:** pairwise cosine distance between A and B, divided by its
  maximum with the implemented numerical guard. The default expression feature
  mask excludes LR genes. The code's documented all-gene fallback applies when
  no non-LR feature remains or the relevant gene metadata are absent.
* **LR compatibility:** form the directed sum
  `S_raw[i,j] = sum(w[L,R] * expression_A[i,L] * expression_B[j,R])` over available
  pairs, then min–max normalize it. A constant sum becomes a zero matrix.
  For nonconstant `S`, the LR cost is `-log(S + 1e-4)`, normalized as in the
  implementation. A constant LR channel contributes zero cost.
* **Geometry:** calculate squared Euclidean distances within each section and
  normalize each distance matrix by its maximum with the implemented guard.
  The Gromov–Wasserstein cost contains the current-plan-dependent cross term
  `-2 * D_A @ P @ D_B.T`; it is recomputed and rescaled at every outer update.

The exact implementation is preserved in
[the bound core](provenance/numerical_validation_snapshot/code/lrot_core.py).
The root `lrot_core.py` is retained for explicitly identified legacy entries;
it is not a substitute for the current source-bound execution protocol.

### 4. Solve using the bound numerical protocol

`paper_protocol.binding.bind` verifies the pinned source files and explicitly
binds the certified balanced solver. `solve_locked` checks the expected input
identity before executing the selected configuration. See
[binding.py](paper_protocol/binding.py) and
[the source manifest](paper_protocol/source_binding.json).

The plan is initialized from uniform marginals, `P = a[:,None] * b[None,:]`.
For each outer update, the implementation constructs

```text
C_total = alpha * C_expression + gamma * C_LR
          + (1 - alpha) * beta * C_GW(P)
effective_epsilon = reg_ot + beta * reg_gw
```

The inner balanced problem starts with at most 1,000 scaling iterations. If
needed, the bound implementation applies reduced dual optimization and its
recorded, bounded feasibility correction. It accepts a plan using the residual
and certificate checks, rather than optimizer status alone. The wrapper's
`inner_max_iter=10000` argument does not mean 10,000 initial scaling iterations.

The locked wrapper uses absolute marginal tolerance `1e-8`, relative marginal
tolerance `1e-5` and relative outer-plan-change tolerance `1e-5`. Acceptance can
occur no earlier than the seventh outer update. Only the wrapper's prescribed
outer-cap retry is permitted; an unaccepted result must not enter a figure or
summary. Case-specific parameters, caps and attempt records remain authoritative.

The stored `loss` is the implementation's diagnostic cost summary. It is not
the full entropy-regularized update objective and its curve alone does not prove
global optimality or convergence. This package does not replace the solver with
an approximate implementation or tune inputs to a historical displayed value.

### 5. Map coordinates and, where applicable, stack sections

For the A-to-B barycentric mapping, row-normalize the plan and calculate
`mapped_A = (P / row_sums(P)) @ coordinates_B`, with the numerical guard used by
the corresponding entrypoint. The reverse mapping uses column normalization
and the transposed plan. They have different spot identities and directions.

Check shapes, finite/nonnegative entries, marginal residuals and zero-mass
rows/columns before interpreting a saved plan. Serial stacking applies the
specified pairwise mappings and transformations in their recorded order; a
single pairwise result does not reproduce an entire 3D experiment.

### 6. Evaluate with the selected metric and direction

| Metric | Definition and boundary |
|---|---|
| Soft label-transfer accuracy / ARI | Row-normalized transport transfers B's known layer labels to A; compare the resulting labels with A. Requires annotations and is not correspondence ground truth. |
| Hard transport matching | Select the maximum-mass target in each row and compare labels. Distinct from soft transfer. |
| Aligned-coordinate nearest-neighbour accuracy | Barycentrically map A, find the nearest original B coordinate and compare labels. Distinct from maximum transport matching. |
| Row entropy | Mean entropy of the row-normalized plan. Lower entropy means greater concentration, not automatically greater accuracy. |
| Same-label transport mass | Fraction of transported mass joining equal labels. Distinct from spot-level accuracy. |
| Case-specific scores and summaries | Use the bound score, mask, denominators, sorting order and aggregation; summary means cannot recover per-spot arrays. |

Do not combine these metrics into one ranking or treat annotation-derived
scores as calibrated probabilities. STAligner embedding evaluation is a separate
pipeline; its complete training protocol and dependencies are separately pinned.
Runtime statements must identify preparation/solve/scoring boundaries and the
actual hardware and repetitions; a single local observation is not a general
speed comparison.

### 7. Save, verify and generate figures

Preserve the input file and ordered-array identities, parameters, seeds, source
hashes, environment, complete attempts/diagnostics, saved plan, metrics and logs.
Generate new outputs in a fresh directory outside the package. Keep historical
results intact and distinguish a new execution from retained historical data.

Figure generation has three separate modes: approved-asset export, supported
drawing from saved results, and case-specific numerical recomputation. The
[figure command table](docs/figure_entrypoints.md) and manifests identify the
scope of each figure. Copying an approved PNG is not an experimental rerun.

## Environment and quick start

The recorded environment is Python 3.10.11 on Windows 11 / AMD64. Use
[the tested environment guide](docs/clean_environment.md) and task-specific
[dependency profiles](docs/environments/). The broader
[requirements.txt](requirements.txt) inventory is not a guarantee that every
external method is prepared. No entrypoint installs dependencies automatically.

Run the following from the code directory. Choose output directories that do
not already exist; the examples below export approved images without solving.

```powershell
python tools/verify_manifest.py --strict
python tools/verify_inputs.py --list
python figures/render_current.py --list
python figures/render_supplement.py --list
python figures/render_current.py --mode asset --output-dir "..\..\lrot_main_assets"
python figures/render_supplement.py --mode asset --output-dir "..\..\lrot_supp_assets"
```

To draw the supported F2 panel from retained records, rather than solve again:

```powershell
python figures/render_current.py --figure F2 --mode records --output-dir "..\..\lrot_F2_records"
```

This retains A–E and the separate 123-value trace; it reconstructs the supported
panel F from saved diagnostic histories. Other figures have different boundaries.
Exact raster identity also depends on the pinned drawing environment and fonts.

## Explicit numerical entrypoints

The following are examples for deliberate numerical work, not automatic setup.
They require their stated data and runtime. This packaging operation executes
none of them.

```text
python experiments/supp/synthetic_seed_robustness.py --verify-only
python experiments/supp/synthetic_seed_robustness.py --run --solver revalidated --output-dir "<new-30-pair-output>"
python experiments/validate_f5_protocol.py --repo-root . --data-dir "<DLPFC-directory>" --output-dir "<new-F5-output>"
python experiments/run_real_visium.py --solver revalidated --input-h5 "<mouse-H5>" --output-dir "<new-mouse-output>"
```

The 30-pair command performs 60 fresh solves when used without checkpoint reuse.
F3/STAligner and F7 have separate preparation, execution and evaluation steps;
follow [the running guide](RUN_REPRO.md) and
[the F7 guide](provenance/f7_prospective/README.md). Do not run a historical
driver's `main` as a shortcut to the bound protocol.

## Resources, attachments and verification scope

Supplementary File S2 is **one CSV with 30 seed pairs**, not 30 CSV files:
[S2](lrot_output/s23_support/Supplementary_File_S2_synthetic_seed_pairs.csv).
Other experiment seeds are recorded in their own configurations and evidence.

The fixed [GO resources](data/go/README.md) are included with their source
records and checksums. Tissue data and CellChatDB binaries remain external.
Retain the resource-specific licence and attribution; project MIT does not
relicense third-party materials. See [external resources](docs/external_resources.md).

Package integrity, approved-image identity, saved-result plotting, independent
saved-array scoring and raw-input-to-solver reproduction are different levels.
Their existing case-specific evidence and limitations remain available in
[the verification guide](docs/release_readiness.md). Not every paper experiment
is certified end to end, and technical repetition is not biological replication.
Historical case identifiers are retained for traceability.

## Licence and snapshot identification

LROT source code is covered by [LICENSE](LICENSE). Identify this local snapshot
using [LOCAL_VERSION.json](LOCAL_VERSION.json), [RELEASE_FILES.txt](RELEASE_FILES.txt)
and the supplied checksum manifests. The download archive has a separate checksum.
[CITATION.cff](CITATION.cff) retains the existing collective contributor entry;
this package does not invent personal authorship, a public hosting address,
a formal software release or an archive DOI.
