# Selected-result verification

This check supplements the earlier [verification scope](release_readiness.md). It changes no computational source, formula, adopted image or experimental result. [Summary](../provenance/final_verification_20261016/summary.json) and [detailed checks](../provenance/final_verification_20261016/checks.json) record the measured scope.

| Item | Check and result |
|---|---|
| Adopted images | All 31 selected images match the corresponding current DOCX media. This is image identity only; final documents remain subject to author revision. |
| Table 1 | Recomputed means, paired differences and Wilcoxon tests from retained 20-seed full-precision metrics; all 24 numeric cells agree. |
| Table 2 | Reconstructed four raw DLPFC pairs and matched source/target input hashes. Eight saved plans reproduce per-spot scores and correctness arrays exactly; displayed correlations and residuals agree. |
| Tables S18-S20 | Fifty raw-input-hash-matched saved plans independently rescored for entropy, residuals, same-layer mass, ARI and accuracy. Displayed values agree, allowing recorded machine-precision residual differences. |
| Tables S2/S8 | Sixty numeric cells agree with the corresponding main-alignment records verified above. |
| Tables S1/S5/S6/S13/S15/S23 | 191 explicitly enumerated numeric cells agree with run outputs or independently aggregated seed/trial statistics. Columns outside the detailed checks are not certified by this count. |
| F5 | Correctness, score-decile accuracies, retention accuracies and threshold coverage/accuracy agree between retained plotting and independent-solve arrays. LROT pooled rho differs at full precision but both round to 0.363. Top-20% accuracy remains LROT 0.7825 / FGW 0.8275. Score and entropy arrays are not identical. |
| F7 | The later same-run payload exactly matches the retained accepted plans. Independent entropy is FGW 5.849532811322741 / LROT 5.817596024603023; reduction 0.5459715801217444%. These support the displayed 5.8495 / 5.8176 / 0.55%. |

F7 panel sources remain distinct: A-D use coordinates, loss history and P from the retained historical plotting payload; E displays manually specified prior coefficients; F reports the later matched-mask summary. Historical P differs from the later accepted P by at most 6.10e-9 per element. Small differences do not establish historical same-run identity. F5/F7 figure captions and Note S25 should retain their applicable provenance boundaries when authors finalize the manuscript.

No new solver or training run was launched for this check. Independent saved-plan scoring, raw-input-identity matching, image identity and fresh optimisation are separate evidence levels. The detailed checks are audit outputs, not substitutes for the underlying full run collection. That large private collection, manuscript text and machine-specific paths are excluded from this release.

Three byte-identical file pairs under `provenance/f3_staligner/` are retained intentionally: embeddings, independently recorded scores and loss histories belong to two complete technical executions. Deleting one would remove the evidence of deterministic repetition. Historical manifests and source-bound numerical files also preserve distinct provenance or executable roles.

Author names/order, public version and repository/current archive identifiers still need confirmation. This candidate retains honest limits and does not claim whole-paper end-to-end reproduction or final submission approval.

## Package code/naming/privacy review

The final archive was checked for missing allowlist members, active figure/source bindings, local project-module names, syntax and local documentation links. Comments/docstrings and release text were scanned for personal names, user/home paths, credentials and conversational remnants. PNG and workbook metadata and 29 pickle/string payloads were also inspected without unpickling or executing them. No user-specific private information or credential was identified by these checks. The generic `D:\LROT_code` documentation example and the upstream STAligner author's public attribution are retained; neither is a private workstation location or user credential. Public authorship/licence information must not be removed as if it were a secret.

One unreferenced historical CellChat log explicitly labelled superseded (2,094-pair convention) is omitted from this delivery; the original remains in local history. The retained selected CellChat summary and all computational sources, arrays and approved images are preserved. Actual archive member paths now use ASCII names, with no case-insensitive collisions or Windows-reserved components. Historical identifiers and provenance dates remain where they describe scientific records. Descriptive historical script numbers are mapped through the figure entrypoint table.

These are scoped code/package checks, not a proof that every external dataset has been bundled or every experimental condition rerun. Raw tissue data, CellChat binaries and the complete run archive remain separately prepared inputs as documented. Formal author/version/repository citation metadata still needs author confirmation.
