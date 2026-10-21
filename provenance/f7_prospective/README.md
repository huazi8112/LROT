# Figure 7 fixed-input reproduction

This is a **new fixed-input reproduction**, not a recovery of the historical
run that produced the approved paper raster. The approved DOCX image remains
unchanged. This protocol used the archived numerical solver source bound by
`paper_protocol/source_binding.json`; the repository's legacy default solver
is a separate implementation.

Place the four original Section 1/2 H5 and positions files in a local data
directory in the relative layout recorded by `reference_input_identity.json`.
They are not distributed in the code package. Their SHA-256 must match that
reference before any calculation. Use the pinned Python environment in
`requirements.txt`; the observed run used Python 3.10.11. Do not adjust the
fixed seeds, barcode order, transform, solver settings, or acceptance threshold
to fit paper numbers.

From the repository root, create one **existing empty** preparation directory; choose fresh nonexistent solver and plot-input directories, and a new image output path; then run:

```text
python experiments/prospective_f7_prepare.py --data-dir <raw-data-dir> --output-dir <prep-dir>
python experiments/prospective_f7_solve_locked.py --prepared-dir <prep-dir> --output-dir <solver-dir>
python figures/build_f7_same_run_plot_input.py --prepared-dir <prep-dir> --solver-dir <solver-dir> --output-dir <plot-input-dir>
python experiments/_render_breast_realcoord_official.py --input-npz <plot-input-dir>/F7_same_run_plot_input.npz --metrics-json <plot-input-dir>/F7_same_run_metrics.json --output-png <new-output-dir>/F7_same_run.png
```

The preparation step checks original file SHA-256, ordered barcodes, ordered
array hashes, LR matrix and source code. The solver runner checks the bound
source fingerprint, input identity, fixed parameters and accepted diagnostics
before reporting scores. The plot-input builder reads only the two accepted
saved plans from **that same new run** and independently computes row entropy.
The final step draws an image from saved results; it does not solve again.

The recorded run used 3,000 ordered spots per section, seeds A=42/B=99,
Section 2 rotation 2 degrees plus translation 0.5, alpha=0.5, beta=0.3,
gamma=0 (FGW) or 0.1 (LROT), and an outer cap of 200. Both methods were
accepted on their first attempt after nine outer updates. Their independently
computed entropies were 5.8495328113 and 5.8175960246 (0.54597158% lower
for LROT), within the predeclared absolute tolerance 0.001 of the retained
summary. The new image is visually equivalent but not byte-identical to the
approved image. Neither image nor paper text was replaced.

`result_manifest.json` identifies local evidence by SHA-256, including the
original failed import log and the successful rerun log. The large prepared
input and saved plans are retained separately and are not part of the release
allowlist. This one case does not certify the whole repository or prove which
historical run produced the approved image. Third-party raw-data redistribution
requires separate permission.
