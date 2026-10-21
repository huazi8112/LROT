# S2D layout repair

The prior D panel was an immutable 1344×1007 raster inserted with an aspect-preserving contain fit into the 1490×1247 cell. This left a smaller plotting region while its text retained the raster's larger relative size. The renderer now plots the six accepted saved diagnostic rows directly into a 1490×1247 panel. D uses the same axes box as E ([180, 97, 1467, 1097]) and the existing S2 font settings; no nonuniform stretching is used.

The plotted accepted outer-iteration values remain 9/39, 9/121 and 9/40 (Full, Shuffled and Adhesion; γ=.1 and γ=1). These are stopping diagnostics, not elapsed times. Earlier S2/D raster exports are excluded from the current release allowlist. The package contains the selected whole-figure asset [assets/S2.png](assets/S2.png), not a second historical S2/D image set. No solver, experiment, parameter or scientific conclusion is changed by this documentation correction.

Run from the repository root:

```text
python figures/render_supplement.py --figure S2 --mode records --output-dir "<new-directory-outside-project>"
```

Released inputs include [the selected S2 asset](assets/S2.png), [S2 saved values](inputs/S2.json), [six accepted D-panel diagnostic rows](inputs/S2D_accepted_diagnostics.json), [the rendering configuration](selected_panel_config.json), and the rendering code identified in [figure_manifest.json](figure_manifest.json). The command writes the regenerated `S2.png`, `S2_D.png`, other panel outputs and validation JSON into the specified new output directory. These generated files are not automatically part of `RELEASE_FILES.txt`; in particular, no standalone D-panel PNG is bundled. Rebuild and verify `MANIFEST.sha256` and `MANIFEST_size.txt` after release-file changes.

The historical D-layout checks are preserved in [S2_layout_repair.json](validation/S2_layout_repair.json) and [S2_saved_record_entry.json](validation/S2_saved_record_entry.json). They identify the pre-title-revision image and retain the paths/hashes recorded at that stage; those historical paths do not imply that the earlier raster files are still distributed. The subsequent title-only revision is recorded in [title_layout.json](validation/title_layout.json). Current selected-image identity is defined by [figure_manifest.json](figure_manifest.json). These records concern saved-result drawing and image identity, not end-to-end experiment reproduction or certification of the author's current working PDF. No drawing or numerical execution is performed by this documentation correction.
