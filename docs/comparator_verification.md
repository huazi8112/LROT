# Targeted comparator and GO verification

This update closes the three selected input/scoring/resource gaps. It does not certify every archived experiment or replace the existing whole-paper verification limits.

| Evidence | Inputs and protocol | Result and limit |
|---|---|---|
| Historical DLPFC STAligner embeddings | 151507→151508 and 151508→151509; 1,000 annotated rows; seeds 42/99; reconstructed barcode/layer order | Independent accuracy .580/.518 and ARI .383938/.315308. Incoming historical training matrices were not recorded. |
| Fresh DLPFC STAligner, selected for Tables S9/S22 | Same prescribed pairs/seeds; normalize_total 10⁴, log1p, PCA50; within-section KNN6 plus reverse edges; MNN50; seed 666; one CPU thread; 500+500 training steps | Accuracy .586/.525, ARI .377807/.302971. Actual training inputs, ordered barcodes, logs, weights and embeddings are bound. One prescribed run per pair; no seed search or donor replication. |
| Selected PASTE/PASTE2 comparisons | 106 retained plans across 23 reconstructed input pairs; full method-specific features/costs/parameters in [experiments.json](../provenance/comparator_audit/experiments.json) and [datasets.json](../provenance/comparator_audit/datasets.json) | 70 directly captured input identities match. The other 36 have input-bound repeat counterparts; largest plan-entry difference is 4.44×10⁻¹⁶, with identical class predictions. Historical capture and later reconstruction remain distinct. |
| Twenty-seed PASTE2 summaries | Seeds 1–20; 1,000+1,000 spots, 150 genes; rotation 5°, batch .15; α=.1; uniform marginals and complete overlap | KL mean .7006, sample SD .073528; Euclidean .4186, SD .141906. Reported rounded summaries remain unchanged. |
| GO BP enrichment | Fixed OBO/GAF; 109 LR genes; 27,431 annotated background genes; foreground-complement Fisher; BH across all 14,584 supported terms | All term records and Figure S3 display rows regenerate exactly. Official current files byte-match the local inputs; original acquisition dates remain unknown. |

Comparator input generation uses the archived two-thread BLAS setting. A seed-8 coordinate reconstruction under one thread differed in floating-point last bits; using the recorded two-thread setting restored exact recorded hashes. This is recorded as an input-generation detail, not hidden by weakening the checksum test. Fresh neural training uses its separately specified one-thread protocol, consistent with the current F3 execution.

## Check retained results

```text
python experiments/verify_comparators.py --output-dir "<new-check>"
```

This independently scores all 106 saved plans, both historical and fresh DLPFC embeddings, and all GO p values/FDR values. It does not run an optimiser or train a model.

## Regenerate comparator inputs or rerun one method

Choose a record ID from `experiments.json`. For example:

```text
python experiments/reproduce_external_comparator.py --record "run_paste2_multiseed_euclid/PASTE2_16.json" --mode inputs --output-dir "<new-synthetic-input-check>"
python experiments/reproduce_external_comparator.py --record "run_dlpfc_baselines/PASTE_1.json" --data-root "<data-root>" --mode inputs --output-dir "<new-DLPFC-input-check>"
```

`--mode run` additionally executes that recorded PASTE/PASTE2 configuration. It requires the external implementation file hash recorded in [environment.json](../provenance/comparator_audit/environment.json); a version label alone is insufficient. PASTE distribution metadata remain unavailable and are not invented. The source is not silently replaced with a newer package. Outputs and differences are retained even when a regression check fails. New timings are not adopted as a speed benchmark.

For fresh neural preparation/training, follow the [DLPFC case guide](../provenance/dlpfc_comparators/README.md). The published preparation entry reproduces all eight actual training-input/order files byte-for-byte. Its training body is the same instrumented pinned implementation used for the two completed fresh runs; only explicit paths and automatic process-local thread settings were added.

## Prepare fixed GO resources

The two fixed ontology and annotation files are bundled in `data/go/` at the repository root, with attribution, source records and checksums. No outer directory or separate GO archive is required. Original acquisition history is not reconstructed from the resource headers.

```text
python tools/prepare_go_resources.py --resource-dir "<unpacked-resource-directory>" --output-dir "<new-resource-copy>"
python tools/prepare_go_resources.py --download --proxy "<your-proxy-if-needed>" --output-dir "<new-download>"
python experiments/go_enrichment.py --obo "<go-basic.obo>" --gaf "<mgi.gaf.gz>" --output-dir "<new-GO-results>"
```

Use one of the first two commands. Downloading is explicit and fails if current upstream bytes change. The GO/MGI official CC BY 4.0 declarations and exact source URLs are retained in [resource_identity.json](../provenance/go_enrichment/resource_identity.json). The current result table is [all_terms.csv](../provenance/go_enrichment/all_terms.csv). The old `lrot_output/lrot_go_enrichment_results.txt` remains a historical result and is not the current paper basis.

The LROT core, numerical source binding, manual priors, all 31 selected figure images and existing paper equations are unchanged.
