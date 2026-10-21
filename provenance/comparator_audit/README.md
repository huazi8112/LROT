# Selected external-comparator input audit

`datasets.json` records ordered expression/spatial/label identities and method-input shapes for 21 simulated seed pairs and two measured DLPFC pairs. `experiments.json` maps each of 106 method records to its dataset, parameters, scores and evidence level. `plans.npz` deduplicates 90 numerical plans; `ordered_labels.npz` supplies the exact scoring order. Source hashes and the retained external environment are recorded separately.

Seventy records directly captured incoming method-input hashes, all verified. The other 36 have numerically equivalent plans in verified input-bound repeats; their historical incoming-input capture remains absent. The largest entrywise difference is 4.44e-16 and class predictions are identical. These statements are deliberately separate. LROT/soft-FGW are not replaced with installed PASTE/PASTE2 implementations.

Input regeneration requires the archived two-thread BLAS setting for exact coordinate hashes. The current independent scoring and input-preparation entries are documented in [the guide](../../docs/comparator_verification.md). Saved-plan verification is not a new optimisation run; reconstructed inputs do not fabricate missing acquisition logs or software distribution metadata.
