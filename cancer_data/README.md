# Breast and mouse raw inputs

Inputs are supplied separately; see the [data preparation guide](../data/README.md) for source identifiers, paths and access links. The historical directory name also covers mouse data and is retained to preserve script paths.

Breast requires two H5 matrices and matching tissue-position CSVs for Block A sections 1/2. Their frozen identities are in `BREAST_SHA256SUMS.txt` and `../data/input_files.json`. Use `python tools/verify_breast_data.py --data-dir "<breast-directory>"` or the unified input checker.

Mouse is CytAssist Fresh Frozen Mouse Brain Post Xenium Rep1 (spatial-exp 2.1.0). Supply its H5 as `mouse_brain.h5`; the frozen SHA-256 is in the input manifest. Run `experiments/run_real_visium.py` with explicit `--input-h5` and a new external `--output-dir`. Revised versus legacy protocols and the imposed-perturbation scoring boundary are described in [RUN_REPRO.md](../RUN_REPRO.md).

These raw files are retained locally where available and excluded from the public release; no data are deleted by packaging.
