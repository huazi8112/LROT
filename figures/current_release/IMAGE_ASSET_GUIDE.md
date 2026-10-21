# Current image assets

`FINAL_PAPER_IMAGE_INDEX.csv` selects exactly 31 manuscript images: F1-F8 and S1-S23. Its paths and hashes identify the approved PNG bytes. This release contains 34 PNG files in total; the other three are the required F8 composition panels in `figures/current_release/components/`.

All eight main figures are in `assets/`. All 23 supplementary figures are in `figures/supplement/assets/`. There is one approved whole image per figure, selected by the CSV.

Superseded whole-figure templates, previous panel exports, redundant F2 diagnostic rasters and unused standalone workflow plots were archived outside the repository. Their numerical records, arrays, saved plans and historical verification records remain preserved. The subsequent layout revision changes path references in selected entrypoints, not computational definitions or experimental values. Historical records can name archived PNGs and old source hashes; those identities describe earlier checks and are not current runtime dependencies.

Current F1 label rendering and supported supplementary local/saved-record drawing use the approved image as their retained-pixel canvas. Only declared regions are replaced; protected pixels are preserved. This deliberately depends on an approved raster and is not full figure reconstruction from raw data. `panels` mode crops the current approved figure on demand; it does not use an older panel image or independently regenerate its numbers.

S23's old-point displacement repair is retired. It was a one-time correction requiring a superseded raster; applying it to the current figure would shift the point twice. The current S23 asset, 30-pair numerical records and original correction evidence remain available. The old stacked supplementary builder now exports/verifies current approved assets only; use the documented records entry for supported numerical drawing.

`RELEASE_FILES.txt` excludes the archived rasters. Legacy experiment scripts may produce standalone diagnostic images when explicitly executed; such new outputs are not automatically selected for the paper or included in a release. Keep new run outputs outside the repository as specified in `RUN_REPRO.md`.

Image identity, saved-result drawing, independent scoring and raw-input execution have different scopes. See `docs/release_readiness.md`; this cleanup does not certify full-repository end-to-end reproduction.

The directory consolidation moves current inputs/assets only. Historical checks retain their original identities; see `provenance/path_relocations.json`.
