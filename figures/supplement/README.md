# Supplementary figure sources

The current manifest selects S1-S23, one approved whole-figure PNG per figure. `asset` exports approved bytes; `records` draws supported regions from saved inputs on the current approved canvas. Outputs must go to a new directory outside the repository.

```text
python figures/render_supplement.py --list
python figures/render_supplement.py --mode asset --output-dir "<new-asset-directory>"
python figures/render_supplement.py --figure S4 --mode records --output-dir "<new-record-drawing-directory>"
python figures/render_supplement.py --figure S2 --mode panels --panel D --output-dir "<new-panel-directory>"
```

| Location | Role |
|---|---|
| `assets/` | Contains all 23 current approved images; no alternative full figures. |
| `inputs/` | Saved records with distinct input/precision/source bindings. |
| `s19_posthoc/` | Saved plan, ordered coordinates/labels, bindings and drawing configuration. |
| `history/`, `validation/` | Historical indices and checks; not current-image selections. |

Old `templates/`, `templates_selected/`, `confirmed_panels/` and the old S19 raster were archived outside the code repository. Supported renderers now retain pixels from the current approved asset. `panels` crops current assets on demand; these crops have current titles and do not represent independent numerical verification. The manifest preserves previous panel-export identities as historical records.

Check `record_render_supported` before records mode. S2/S4/S5/S7/S8/S9/S15/S17/S19 retain saved-record drawing support. Patch mode supports S3/S6/S12/S20/S21 with the exact pinned Arial fonts. S23's one-time old-point correction is retired; use asset export for its current figure. Numerical records and correction evidence remain preserved. Recorded drawing environment: Python 3.10.11, Matplotlib 3.10.7, NumPy 2.2.6, Pillow 12.0.0. A records byte mismatch fails and never authorizes replacing a manuscript image.

S19B recomputes mapping from the saved plan and ordered coordinates without loading H5, expression or a solver. S4A uses retained four-decimal summaries and protects B/C. S2D remains bound to six accepted saved diagnostics. S7 separates outer solve time from LR preparation. S8 retains both tested PASTE2 configurations. S15 C/D use sample SD over five unique seeds; B shows saved means without added error bars.

See [image roles](../current_release/IMAGE_ASSET_GUIDE.md), [running guide](../../RUN_REPRO.md) and [verification scope](../../docs/release_readiness.md). No full-repository end-to-end certification is claimed.
