# Verification scope

Package and image identities are recorded in `LOCAL_VERSION.json` and the figure manifests. Package integrity and approved-image export do not certify numerical reproduction. Document/PDF hashes retained in the evidence identify their stated historical baselines; manuscripts are not included and this local packaging operation does not certify a later manuscript.

| Case | Evidence and boundary |
|---|---|
| F1 | Approved schematic asset retained; full schematic redrawing is a separate output. |
| F2 | Four 80-update diagnostic curves independently reproduced. The 123-value main trace is retained, not newly executed. No convergence or historical recovery claim. |
| F3 | All 30 controlled pairs freshly solved (60 accepted plans) and rescored to the current CSV/F3 values. The adopted STAligner fixed case is .930, ARI .853503; complete deterministic technical repeats match. Historical .934 is superseded, not recovered. |
| F4/F8 | Saved-result drawing available; this does not certify every raw-data condition. |
| F5 | Fresh Br5292 raw-input case with four accepted solves and per-spot scoring. Historical per-spot run binding is not recovered; the redraw differs from the approved raster. Donor-2 certification is not implied. |
| F6/mouse | Fixed raw-input perturbation case independently executed. Approved raster and saved-result drawing have separate identities. Region accuracy is unavailable. |
| F7 | Fresh fixed-input preparation and two solves plus independent scoring are recorded. The approved raster uses an older saved plan and later matched-mask summary; its historical same-run origin is not authenticated. |
| S1–S23 | All approved assets exportable; saved-record/reconstruction support is figure-specific. S19B is a post-hoc saved-plan reconstruction. All supplementary conditions are not certified end to end. |

Numerical evidence: [fresh acceptance](../provenance/reproduction_acceptance/summary.json). Entry commands: [figure entrypoint table](figure_entrypoints.md) and [RUN_REPRO.md](../RUN_REPRO.md). Resource interpretation: [lr_resources.md](lr_resources.md).

Using this package:

1. Distribute the allowlist plus its two manifests; verify the clean directory strictly. Raw inputs, private backups, manuscripts and old Git history remain excluded.
2. Prepare hash-matched external inputs as described in [data/README.md](../data/README.md). CellChatDB now has a verified matching public commit (declared GPL-3) and byte-exact PKL reconstruction, documented in [fixed external resources](external_resources.md). Its binaries remain outside the MIT code package. Original acquisition/freeze dates and other unresolved resource records remain explicitly unknown.
3. Keep the pinned STAligner dependency requirements and technical-repeat scope. The bundled upstream MIT source plus disclosed patch reconstruct the exact six training files offline; follow [fixed external resources](external_resources.md). No automatic upstream substitution or backend change is permitted.
4. When comparing a manuscript with the included assets, identify its exact DOCX/PDF files and check figure/caption consistency. File-package hashes alone do not establish document synchronization.
5. Identify this local snapshot by the supplied file manifests. `CITATION.cff` retains its existing collective contributor entry; no public hosting address or archive DOI is assigned here.

No full-repository end-to-end certification or general comparator ranking is claimed. Local preparation does not perform or authorize a remote push.

Fixed external-source resolution: matching public CellChatDB RDA acquisition and exact derived PKL rebuilding passed. Official STAligner source plus the disclosed patch reconstructed all six training files; one complete source-isolated 1000-step run reproduced .930, ARI .853502856, all losses, returned embedding and final weight tensors. The recorded dependency environment was reused; no fresh dependency installation is certified. See [source reconstruction evidence](../provenance/f3_staligner/source_rebuild_verification.json).

Subsequent clean-dependency verification: a new venv with no system/user site packages installed the fixed dependencies from public package sources, then completed the fixed 1000-step F3 run. Accuracy .930, ARI .853502856, returned embedding, all losses and final weight tensors match. Earlier statements about dependency reuse describe the earlier check only. See [clean-environment evidence](../provenance/clean_environment/verification.json). Cross-platform and whole-paper certification remain outside this scope.

Selected-result update: [Table 1/Table 2, S18-S20 and F5/F7 checks](final_verification.md) extend the saved-plan verification evidence. Table 2 now includes raw-input-bound scoring for both donors; this does not imply new donor-2 optimisation or all-condition certification.

## Independent repository distribution

The current allowlist includes frozen GO/MGI inputs and the S1 workbook inside this repository. No sibling resource/attachment directory is needed. Historical companion-package and earlier inventory records remain unchanged evidence of their own scope; `RELEASE_FILES.txt` and the current manifests define this package. Current official DLPFC donor names are mapped to retained run aliases in [data/dlpfc_sample_aliases.json](../data/dlpfc_sample_aliases.json). This packaging/naming change does not revalidate all experiments or resolve the biological-prior, authorship or formal archival-release questions.
