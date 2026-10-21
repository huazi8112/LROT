# Fixed external resources

The machine-readable source and checksum record is [external_sources.json](../data/external_sources.json). Original experiment inputs, parameters and results are unchanged. The public commit selected here is a verified reconstruction source, not proof of the original download or installation date.

## CellChatDB.mouse

The RDA at the fixed [sqjin/CellChat commit e4f68625](https://github.com/sqjin/CellChat/tree/e4f68625b074247d619c2e488d33970cc531e17c) is byte-identical to the retained RDA. That commit's DESCRIPTION declares CellChat 1.6.1 and GPL-3; the complete [upstream licence](../data/cellchat_LICENSE.txt) is retained. The newer jinworks database checked during this review differs and is not substituted.

The code package excludes both RDA and PKL. Obtain the fixed public RDA separately, then rebuild PKL using the hash-bound legacy conversion rules. The preparation tool preserves symbol handling, receptor-component projections and .9/.6 weights. Its fixed 150-gene list was regenerated using the original builder protocol; it is distinct from the 200-gene synthetic figure archive. Reconstruction produced the same 2059 full pairs, 38 selected pairs and exact PKL bytes. Do not interpret these implementation constants as database confidence.

Explicit network preparation (no manuscript or local file content is transmitted):

```powershell
python tools/prepare_external_resources.py cellchat --download --output-dir ../lrot_external/cellchat
python tools/verify_inputs.py --dataset cellchat --resource-dir ../lrot_external/cellchat
```

If a proxy is required, add `--proxy "<proxy-url>"` to the preparation command using your own network configuration. The option applies only to that request and does not change system settings. For offline preparation, replace `--download` with `--rda "<path-to-pinned-CellChatDB.mouse.rda>"`. Dependencies include the recorded `rdata==0.11.2` and its dependencies. No dependency installer is invoked. Output must be a new external directory; originals are not overwritten.

The preparation record identifies the public commit and generated checksums. The unknown historical acquisition date remains unknown. GPL-3 is the upstream declaration; this resource and its derived PKL are not relicensed by LROT's MIT licence or included as MIT data. Separate acquisition avoids adding these binaries to the code distribution.

## STAligner fixed F3 source

Five of the six recorded training files are identical to [zhoux85/STAligner commit 85ca4ce3](https://github.com/zhoux85/STAligner/tree/85ca4ce34072a37b949733aa55a3d3854ed89038). The remaining file, `mnn_utils.py`, differs by the disclosed [patch](../provenance/f3_staligner/source_recipe/mnn_utils.patch): optional annoy/hnswlib imports, fallback to the existing exact sklearn neighbour function, and the explicit `n_neighbors` keyword. The full [MIT licence and upstream copyright notice](../provenance/f3_staligner/source_recipe/LICENSE) accompany the six upstream files. The original local patch author/date and installation route were not recovered; no authorship claim is made.

Reconstruct the exact six training files offline:

```powershell
python tools/prepare_external_resources.py staligner --output-dir ../lrot_external/staligner
```

This checks upstream hashes, applies the disclosed patch and checks all six resulting hashes against the recorded environment. It also generates metadata for this newly reconstructed source tree, so the original STAligner installation is not needed to supply package metadata. It does not install Python or other dependencies.

With the [recorded dependency environment](../provenance/f3_staligner/environment.json) already available, use the reconstructed source:

```powershell
$env:PYTHONPATH = (Resolve-Path ../lrot_external/staligner).Path
python experiments/reproduce_f3_staligner.py --mode train --output-dir ../lrot_f3_source_check
```

The original entrypoint enforces all six source hashes, package versions, single-thread execution and the exact-neighbour backend. Approximate annoy/hnswlib backends must be absent for this fixed case; the gate rejects them. No gate is bypassed, seed changed or training protocol altered. The command requires a new output directory and refuses network connections during training. It checks the complete 1000-step execution against the fixed saved embedding and scores.

Source reconstruction and reuse of recorded dependencies are different from a new clean dependency installation. This review does not certify an independent clean installation or cross-platform equality. See [F3 verification scope](../provenance/f3_staligner/README.md). Historical .934 is not recovered, and none of this certifies every experiment in the paper.

Subsequent clean-dependency verification: a new venv with no system/user site packages installed the fixed dependencies from public package sources, then completed the fixed 1000-step F3 run. Accuracy .930, ARI .853502856, returned embedding, all losses and final weight tensors match. Earlier statements about dependency reuse describe the earlier check only. See [clean-environment evidence](../provenance/clean_environment/verification.json). Cross-platform and whole-paper certification remain outside this scope.
