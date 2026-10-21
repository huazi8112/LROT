# Code layout

For figure IDs, commands and input/scope mappings, use the [figure entrypoint table](figure_entrypoints.md) ([中文](figure_entrypoints_zh.md)).

| Directory | Purpose |
|---|---|
| `data/` | Unified input checksums and LR/GO resource instructions; not an additional tissue dataset. |
| `cancer_data/`, `dlpfc_data/` | External raw-input instructions and retained compatibility paths. The historical cancer path also covers mouse brain. |
| `experiments/`, `experiments/supp/` | Experiment/preparation/verification entrypoints and supplementary cases. Some historical helpers draw figures; use the documented primary commands. |
| `figures/current_release/assets/` | Exactly eight current main figures, F1-F8. |
| `figures/current_release/inputs/` | Saved main-figure payloads, array/binding files and shared fixed-case inputs. |
| `figures/current_release/components/` | Three required F8 composition panels. |
| `figures/current_release/evidence/` | Retained F2 composition evidence. |
| `figures/supplement/assets/` | Exactly 23 current supplementary figures, S1-S23. |
| `figures/supplement/inputs/` | Saved supplementary values and bindings. |
| `figures/supplement/s19_posthoc/` | Saved plan and ordered inputs for explicitly post-hoc S19 mapping. |
| `figures/supplement/history/`, `validation/` | Historical manifest and drawing checks; not alternative current figures. |
| `lrot_output/` | Retained numerical summaries; see its README for six clearly named result subdirectories. |
| `paper_protocol/` | Current fixed numerical execution API and source-identity checks. |
| `provenance/` | Case inputs, environments, scores, source identities and acceptance records. |
| `tools/` | Read-only checks, resource preparation and release construction. |

The numerical snapshot in `provenance/numerical_validation_snapshot/code/` is an active source dependency of `paper_protocol.binding`, despite its archival name. `archived_wrapper.py` is likewise imported by the current protocol; neither is an unused backup. The four computational snapshot sources and extracted numerical definitions retain their original hashes. They are not interchangeable with the legacy root solver.

The raw-input directory names and legacy script filenames are retained where fixed source/input bindings depend on them. In particular, some plotting filenames use older manuscript figure numbers; the manifests and 31-image index define current F/S identities. They must not be inferred from those filenames.

Current top-level figure commands are `figures/render_current.py` and `figures/render_supplement.py`. Asset export and supported retained-record drawing require new output locations outside the repository. Historical sources/checks preserve the paths and hashes that were actually recorded; [path_relocations.json](../provenance/path_relocations.json) maps relocated files without rewriting evidence.

See [RUN_REPRO.md](../RUN_REPRO.md), [input preparation](../data/README.md), [package contents](package_contents.md) and [verification scope](release_readiness.md).
