# Figure entrypoints

[中文](figure_entrypoints_zh.md)

Run from the repository root using the recorded environment. This table lists supported primary entries, not every legacy helper. Replace `<out>` with a fresh directory outside the repository; use a different directory for each invocation. Replace `<S>` with one figure ID from its row. Font placeholders require local files matching the manifest hashes. Raw-input and dependency placeholders refer to separately prepared, hash-matched resources; no command downloads them automatically.

Check package integrity first, then use the copy-only commands for exact approved images. These commands do not run experiments:

```text
python tools/verify_manifest.py
python figures/render_current.py --mode asset --output-dir "<new-main-assets>"
python figures/render_supplement.py --mode asset --output-dir "<new-supplement-assets>"
```

Inputs below are primary inputs, not exhaustive file lists. Complete executable pins are in the [main manifest](../figures/current_release/figure_manifest.json) and [supplement manifest](../figures/supplement/figure_manifest.json); approved image identities are in the [31-image index](../figures/current_release/FINAL_PAPER_IMAGE_INDEX.csv). Source filenames can contain historical figure numbers; use the figure IDs in this table.

| Figure / task | Command | Primary inputs | Verification scope |
|---|---|---|---|
| F1 | `python figures/render_current.py --figure F1 --mode records --output-dir "<out>"` | [F1.png](../figures/current_release/assets/F1.png) | Replaces the declared LR-prior label band; retains other raster pixels. Schematic, not an experimental rerun. |
| F2 | `python figures/render_current.py --figure F2 --mode records --output-dir "<out>"` | [80-update histories](../figures/current_release/inputs/F2_four_gamma_histories_current_code.npz); [F2.png](../figures/current_release/assets/F2.png) | Draws panel F from saved histories; retains A–E and the 123-value main trace. No historical recovery or convergence certification. |
| F3 | `python figures/render_current.py --figure F3 --mode records --output-dir "<out>"` | [F3.json](../figures/current_release/inputs/F3.json) | Draws saved metrics and paired-case statistics. STAligner scoring/training and 30-pair solving are separate entries below. |
| F4 | `python figures/render_current.py --figure F4 --mode records --output-dir "<out>"` | [DLPFC NPZ](../lrot_output/arrays/lrot_data_dlpfc.npz); [saved summary](../lrot_output/dlpfc_results/lrot_dlpfc_results.txt) | Saved-array/result drawing only; no raw-DLPFC computation. |
| F5 | `python figures/render_current.py --figure F5 --mode records --output-dir "<out>"` | [per-spot NPZ](../figures/current_release/inputs/F5_current_protocol_plot_only.npz); [input binding](../figures/current_release/inputs/F5_input_binding.json) | Draws saved per-spot results. This redraw is not byte-identical to the approved image; use asset export for approved bytes. |
| F6 | `python figures/render_current.py --figure F6 --mode records --output-dir "<out>"` | [mouse NPZ](../lrot_output/arrays/lrot_data_real_visium.npz) | Saved-result drawing. Measured expression with generated coordinates/perturbation; region accuracy is not applicable. |
| F7 | `python figures/render_current.py --figure F7 --mode records --output-dir "<out>"` | [historical plan payload](../figures/current_release/inputs/F7_historical_plan_plot_only.npz); [binding](../figures/current_release/inputs/F7_source_binding.json) | Draws the retained plan plus later matched-mask summary. Historical same-run origin remains unbound; new fixed-input computation is a separate entry. |
| F8 | `python figures/render_current.py --figure F8 --mode records --output-dir "<out>"` | [F8.json](../figures/current_release/inputs/F8.json); [three panels](../figures/current_release/components/) | Composes saved panels; does not rerun 3D alignment or multi-seed experiments. |
| S2, S4, S5, S7, S8, S9, S15, S17, S19 | `python figures/render_supplement.py --figure "<S>" --mode records --output-dir "<out>"` | [supplement manifest](../figures/supplement/figure_manifest.json); [saved inputs](../figures/supplement/inputs/); [S19 plan/bindings](../figures/supplement/s19_posthoc/) | Saved-record drawing on approved canvases. S4A uses rounded summaries; S19B reconstructs from a saved plan. Other retained panels are not independently regenerated. |
| S3, S6, S12, S20, S21 | `python figures/render_supplement.py --figure "<S>" --mode patch --font-regular "<arial.ttf>" --font-bold "<arialbd.ttf>" --output-dir "<out>"` | [saved patch records](../figures/supplement/inputs/); [pinned font hashes](../figures/supplement/figure_manifest.json) | Local result/label revisions using the exact pinned fonts; no new experiment. Not full reconstruction from raw inputs. |
| S1, S10, S11, S13, S14, S16, S18, S22, S23 | `python figures/render_supplement.py --figure "<S>" --mode asset --output-dir "<out>"` | [approved S assets](../figures/supplement/assets/) | Approved-image export through this primary entrypoint only. S23 historical point-displacement patch is retired. |
| F2: new curves | `python figures/reconstruct_f2_curves_current.py --output-dir "<out>"` | [synthetic archive](../lrot_output/arrays/lrot_data_synthetic.npz); [fixed design](../experiments/synthetic_design_manifest.json) | Executes four fixed 80-update diagnostic runs. The separate 123-value main trace remains retained. Not convergence proof or historical-run recovery. |
| F3: STAligner score | `python experiments/reproduce_f3_staligner.py --mode verify --output-dir "<out>"` | [ordered inputs and two saved embeddings](../provenance/f3_staligner/) | Scores retained embeddings without training. Fixed-case agreement .930, ARI .853503; not recovery of historical .934. |
| F3: STAligner train | `python experiments/reproduce_f3_staligner.py --mode train --output-dir "<out>"` | [fixed case and environment pins](../provenance/f3_staligner/README.md); [offline source preparation](external_resources.md); [dependency pins](../provenance/f3_staligner/environment.json) | One complete deterministic 1,000-step fixed-case training run; exact external STAligner sources/dependencies required. No download, seed search or whole-paper certification. |
| F3C / S23: 30 pairs | `python experiments/supp/synthetic_seed_robustness.py --verify-only`<br>`python experiments/supp/synthetic_seed_robustness.py --run --solver revalidated --output-dir "<out>"` | [registered pairs and summaries](../lrot_output/s23_support/); [bound numerical sources](../paper_protocol/source_binding.json) | First command checks retained records without solving. Second performs 60 fresh solves for 30 registered pairs when using a new directory without --resume. Historical records remain separately identified. |
| F5: Br5292 raw case | `python experiments/validate_f5_protocol.py --repo-root . --data-dir "<dlpfc-dir>" --output-dir "<out>"` | [frozen input binding](../figures/current_release/inputs/F5_input_binding.json); [external DLPFC files](../data/README.md) | Four solves plus independent per-spot scoring for the fixed Br5292 case. Does not recover historical run binding or certify donor 2. |
| F6: mouse raw case | `python experiments/run_real_visium.py --solver revalidated --input-h5 "<mouse.h5>" --output-dir "<out>"` | [external mouse H5/checksum](../data/README.md); [input-order evidence](../provenance/mouse_perturbation/) | Recomputes the fixed expression/perturbation case with the recorded input precision and bound solver. Generated coordinates; no region accuracy. |
| F7: new fixed-input case | [Four-step preparation → solve → plot-input → drawing commands](../provenance/f7_prospective/README.md) | [raw-file/order bindings](../provenance/f7_prospective/reference_input_identity.json); [external breast files](../data/README.md) | Two fresh solves after raw-input verification; scoring/plotting use the same new plans. Approved historical raster is separate. Preparation requires an existing empty directory; later stages use fresh destinations. |

Image export, retained-record drawing, saved-embedding scoring and fresh numerical execution have different scopes. A drawing match does not certify every experimental condition. Use the [running guide](../RUN_REPRO.md) for setup, [F7 guide](../provenance/f7_prospective/README.md) for its multi-stage directory rules, and [verification scope](release_readiness.md) / [numerical acceptance evidence](../provenance/reproduction_acceptance/summary.json) for limits. This documentation revision runs no experiment and replaces no paper image.


## Targeted checks

| Paper item | Entry | Scope |
|---|---|---|
| Tables S9/S22 | `experiments/prepare_dlpfc_staligner.py` → `experiments/train_dlpfc_staligner.py` | Two explicitly input-bound fresh STAligner cases; preparation and training use separate processes. |
| External comparisons and GO | `python experiments/verify_comparators.py --output-dir "<new-out>"` | Independent scoring of 106 retained plans, four historical/fresh DLPFC embeddings and 14,584 GO term statistics. No optimiser or training. |
| Figure S3 statistics / Note S16 | `python experiments/go_enrichment.py --obo "<go-basic.obo>" --gaf "<mgi.gaf.gz>" --output-dir "<new-out>"` | Regenerates the complete term table from exact resource checksums; approved image bytes remain unchanged. |

Commands, inputs and evidence levels: [targeted verification guide](comparator_verification.md).
