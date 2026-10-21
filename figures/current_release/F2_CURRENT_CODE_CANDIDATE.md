# Adopted Figure 2 diagnostic curves

The selected final manuscript adopts the recorded current-code four-curve candidate in panel F. Panels A-E and the 123-value main trace are retained. Four ablation curves use 80 outer updates each, with gamma 0, 0.05, 0.1 and 0.2. The selected gamma=.1 curve equals the first 80 values of the retained main trace.

Independent execution reproduced these curves with a maximum absolute difference of 2.22e-16. This is a legacy diagnostic protocol, not proof of convergence, fixed-objective descent, or historical ablation recovery. Figure 3C remains the primary revalidated comparison.

Run `python figures/reconstruct_f2_curves_current.py --output-dir <new-directory-outside-code>` to calculate curves. Run `python figures/render_current.py --figure F2 --mode records --output-dir <new-directory-outside-code>` to redraw the adopted composite from pinned saved histories. Both commands preserve approved files. See `figure_manifest.json` for current hashes; older candidate entries are historical.
