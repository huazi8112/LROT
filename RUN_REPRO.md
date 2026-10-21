# Running and reproduction

Choose a task using the [figure entrypoint table](docs/figure_entrypoints.md), then follow the setup and evidence boundaries below.

Run commands at the repository root. Use Python 3.10.11. For a new installation follow [the tested task profiles and lock](docs/clean_environment.md); `requirements.txt` retains the broader research inventory and now includes the two previously omitted fixed PyG CPU extensions. Existing numerical version pins are unchanged. Full STAligner training additionally requires the external source bytes and versions in `provenance/f3_staligner/environment.json`. No tool below installs or downloads dependencies automatically.

## 1. Prepare and check inputs

[data/README.md](data/README.md) lists dataset identifiers, access pointers, paths and checksum scope. Keep raw tissue inputs and CellChatDB binaries outside the public package. Fixed GO/MGI files are included in `data/go/`. For external tissue inputs, use a data root containing `cancer_data/` and `dlpfc_data/`:

```text
python tools/verify_repository_package.py
python tools/verify_inputs.py --dataset go --data-root .
python tools/verify_inputs.py --list
python tools/verify_inputs.py --dataset breast --dataset dlpfc --data-root "<data-root>"
python tools/verify_inputs.py --dataset cellchat --resource-dir "<authorized-resource-directory>"
python tools/verify_manifest.py
```

Choose only the datasets needed for a run. Missing or changed files fail the input check. Newly inventoried hashes identify available local bytes; they do not recover missing historical acquisition records. CellChatDB resources remain separately supplied; follow the recorded upstream GPL-3 notice and fixed preparation recipe.

## 2. Export or redraw figures

Output directories must be new and outside the repository. Replace placeholders with actual paths.

```text
python figures/render_current.py --list
python figures/render_supplement.py --list
python figures/render_current.py --mode asset --output-dir "<new-main-assets-directory>"
python figures/render_supplement.py --mode asset --output-dir "<new-supplement-assets-directory>"
python figures/render_current.py --figure F3 --mode records --output-dir "<new-F3-drawing-directory>"
python figures/render_supplement.py --figure S4 --mode records --output-dir "<new-S4-drawing-directory>"
```

Asset mode exports all eight main or 23 supplementary approved images without solving. Records mode draws supported figures from pinned saved inputs; it is not a new numerical experiment. Main figures default to records mode. Figure/input hashes are checked before use. See [supplementary modes](figures/supplement/README.md) and the [figure index](lrot_output/FIGURE_INDEX.md).

S4A uses retained four-decimal summaries and preserves B/C pixels. S19B uses a saved transport plan and ordered coordinates/labels; the current approved image supplies protected pixels outside the redrawn B interior. Font/runtime differences can affect drawing bytes. Failure never authorizes replacing the approved manuscript image.

## 3. Reproduce numerical cases

Use `paper_protocol/` for revised numerical execution; explicit legacy modes remain historical. Keep inputs, precision, seeds and solver binding fixed.

```text
python experiments/supp/synthetic_seed_robustness.py --verify-only
python experiments/supp/synthetic_seed_robustness.py --run --solver revalidated --output-dir "<new-30-pair-directory>"
python experiments/reproduce_f3_staligner.py --mode verify --output-dir "<new-F3-score-directory>"
python experiments/reproduce_f3_staligner.py --mode train --output-dir "<new-F3-training-directory>"
python experiments/run_real_visium.py --solver revalidated --input-h5 "<mouse-H5>" --output-dir "<new-mouse-directory>"
python experiments/validate_f5_protocol.py --repo-root . --data-dir "<DLPFC-directory>" --output-dir "<new-F5-directory>"
```

The 30-pair `--run` performs 60 fresh solves without checkpoint reuse. F3 `verify` scores saved embeddings; `train` performs complete deterministic training and rejects changed external source/dependency bindings. F5 performs four accepted solves for the Br5292 case. Mouse uses measured expression, generated grid coordinates and imposed perturbations; region accuracy is not applicable. Its revised float16-to-float32 input protocol is unchanged.

F2 curves: `python figures/reconstruct_f2_curves_current.py --output-dir "<new-F2-directory>"`. Four 80-update diagnostic histories are calculated; the separate 123-value trace is retained, not claimed freshly reproduced. F7 preparation/execution and scoring bindings are documented in [the case guide](provenance/f7_prospective/README.md).

## 4. Interpret checks and package outputs

[docs/release_readiness.md](docs/release_readiness.md) is the single summary of current evidence and limitations; `provenance/reproduction_acceptance/summary.json` preserves fresh numerical acceptance. Do not interpret package integrity, image equality or technical repeats as full-repository certification or biological replication.

`RELEASE_FILES.txt` defines distribution. Build a new package with `python tools/build_release.py --output "<new-package-directory>"`, then run its `tools/verify_manifest.py --strict`. Supplementary File S1 is bundled in `supplementary_files/LROT_Supplementary_Tables_S11_S12.xlsx`; the selected 30-row Supplementary File S2 CSV is also bundled. These attachment identifiers are distinct from Figures S1 and S2. New run outputs and raw resources are not automatically added to the package.

[中文](RUN_REPRO_zh.md)

## Fixed external-source preparation

See [external_resources.md](docs/external_resources.md) for pinned public CellChatDB acquisition/exact PKL rebuilding and offline reconstruction of the six fixed STAligner training files. The original training source/version/backend gates remain enforced. This source recipe does not certify a new clean dependency installation or recover historical .934.


Targeted verification: see [comparator and GO audit](docs/comparator_verification.md). This adds independent scoring/input records and two fresh DLPFC neural runs; it is not whole-paper certification.
