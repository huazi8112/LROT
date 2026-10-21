# Saved numerical results

This directory contains retained numerical summaries and plotting arrays. Current approved images are stored under `figures/current_release/assets/` (F1-F8) and `figures/supplement/assets/` (S1-S23), with three required F8 composition panels under `figures/current_release/components/`. There are no PNG files in this results directory after the layout revision.

| Location | Contents |
|---|---|
| `arrays/` | Five retained NPZ plotting/input archives; not raw tissue downloads. |
| `dlpfc_results/` | Saved DLPFC text read by the current F4 renderer. |
| `dlpfc_donor2/` | Second-donor summaries. |
| `s17_support/` | Distinct ceiling/epsilon scans and LR-program controls; filenames and values preserved. |
| `s23_support/` | Thirty-pair records, registry and the selected supplementary CSV attachment. |
| Root TXT/JSON/CSV files | Retained summaries, diagnostics and machine/human-readable tables. |

The former Chinese intermediate directories and duplicated S17 directory levels were consolidated without deleting a result. Historical records retain original paths and hashes; [path relocations](../provenance/path_relocations.json) resolve their locations. `lrot_paths.py` searches current locations first and retains compatibility with externally supplied pre-migration layouts.

The authoritative image selection is [the 31-image index](../figures/current_release/FINAL_PAPER_IMAGE_INDEX.csv). Earlier image counts and retained-template descriptions are superseded. Old rasters remain in the private project archive, outside this repository. Legacy scripts can produce diagnostic images when explicitly run; outputs are not automatically paper selections or release files.

Use [the running guide](../RUN_REPRO.md) and match input/protocol identities when comparing results. Asset export, saved-result drawing, independent scoring and raw-input execution have different scopes; directory cleanup does not establish full experiment reproduction.

The superseded 2,094-pair CellChat reconstruction log is omitted from this delivery. Its original and earlier relocation record remain in the local historical archive. The retained selected summary is `lrot_cellchat_results.txt`; fixed-resource reconstruction evidence is documented in [external resources](../docs/external_resources.md). No scientific result or source implementation was changed by this packaging cleanup.
