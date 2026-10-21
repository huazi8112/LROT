# External data and resource preparation

Raw tissue inputs and CellChatDB binaries are not bundled. Fixed GO/MGI inputs are bundled in `data/go/`. This guide and `input_files.json` replace old numbered-section pointers. Prepare only the inputs needed for the selected experiment; preserve file names and barcode order.

| Dataset | Required location relative to an external data root | Access/source |
|---|---|---|
| Breast, Block A sections 1/2 | `cancer_data/V1_Breast_Cancer_Block_A_Section_{1,2}.h5` and matching section directories containing `spatial/tissue_positions_list.csv` | [Section 1](https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_web_summary.html), [Section 2](https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_2/V1_Breast_Cancer_Block_A_Section_2_web_summary.html); 10x spatial-exp 1.1.0 |
| Mouse brain | `cancer_data/mouse_brain.h5` | [CytAssist Fresh Frozen Mouse Brain Post Xenium Rep1 H5](https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/CytAssist_FreshFrozen_Mouse_Brain_Post_Xenium_Rep1/CytAssist_FreshFrozen_Mouse_Brain_Post_Xenium_Rep1_filtered_feature_bc_matrix.h5); spatial-exp 2.1.0, not olfactory bulb |
| DLPFC | `dlpfc_data/{sample}_filtered_feature_bc_matrix.h5`, `{sample}_tissue_positions_list.txt`, and `barcode_level_layer_map.tsv` | Samples 151507/151508/151509 and 151669/151670/151671; [spatialLIBD H5 service](https://spatial-dlpfc.s3.us-east-2.amazonaws.com), [HumanPilot coordinates/labels](https://github.com/LieberInstitute/HumanPilot/tree/master/10X) |
| GO/MGI | `data/go/go-basic.obo`, `data/go/mgi.gaf.gz` | Recorded headers: GO releases/2026-07-26; GAF 2.2, generated 2026-05-21. [GO access](https://current.geneontology.org/ontology/go-basic.obo), [MGI access](https://current.geneontology.org/annotations/mgi.gaf.gz) are mutable current links, not guaranteed historical bytes. |
| CellChatDB | An authorized resource directory containing `CellChatDB.mouse.rda` and `CellChatDB_mouse_lr.pkl` | The fixed sqjin/CellChat commit e4f68625 supplies an RDA byte-identical to the retained input; DESCRIPTION declares 1.6.1 and GPL-3. Follow [fixed resource preparation](../docs/external_resources.md) for the pinned URL and exact PKL reconstruction. Do not substitute the newer database. |

DLPFC direct H5 URLs use `https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/{sample}_filtered_feature_bc_matrix.h5`; HumanPilot paths are `10X/{sample}/tissue_positions_list.txt` and `10X/barcode_level_layer_map.tsv`. The existing download helper is an optional access aid; it is not an integrity check. This revision does not run it or change its behavior.

```text
python tools/verify_inputs.py --list
python tools/verify_inputs.py --dataset breast --dataset mouse --dataset dlpfc --dataset go --data-root "<data-root>"
python tools/verify_inputs.py --dataset cellchat --resource-dir "<authorized-resource-directory>"
```

The checker reads files, prints missing/mismatched inputs and returns nonzero on failure. It never downloads, edits inputs or runs a solver. `input_files.json` distinguishes previously bound identities from checksums first inventoried for this packaging revision. A match establishes file identity only, not historical acquisition, licensing or scientific reproduction.

Point numerical entrypoints explicitly to their data: mouse `--input-h5`, F5 `--data-dir`, and F7 options in its case guide. Legacy loaders may use `LROT_DATA_DIR` (breast) or `LROT_DLPFC_DATA_DIR`; portable batch preparation requires `LROT_DATA_ROOT`. Do not assume the input checker sets these variables.

This repository includes the fixed GO/MGI resources and their attribution in `data/go/`; original tissue inputs and the two CellChatDB binaries remain external. Matching CellChatDB public acquisition and exact derived PKL reconstruction are now documented in [external_resources.md](../docs/external_resources.md). This identifies an available matching source, not the historical acquisition date. Other missing historical resources must be reported, not silently replaced. The PKL builder's .9/.6 coefficients are implementation constants; manually specified 61/68-entry priors are separate. See [resource evidence](../docs/lr_resources.md) and the preserved `resource_provenance.json` for biological and source limitations. The fixed GO/MGI files now independently regenerate every reported term and match official current downloads; original acquisition dates remain unknown. See [GO verification](../provenance/go_enrichment/README.md).


Targeted verification: see [comparator and GO audit](../docs/comparator_verification.md). This adds independent scoring/input records and two fresh DLPFC neural runs; it is not whole-paper certification.

## Bundled resource and attachment checks

Run `python tools/verify_inputs.py --dataset go --data-root .` from this repository root to check the included OBO/GAF bytes without downloading or analysing them. The resource headers and original acquisition-history limits are retained in [data/go/README.md](go/README.md). Supplementary File S1 is included [here](../supplementary_files/LROT_Supplementary_Tables_S11_S12.xlsx); Supplementary File S2 remains [the single selected 30-pair CSV](../lrot_output/s23_support/Supplementary_File_S2_synthetic_seed_pairs.csv). For tissue files, supply the external data root and use only the required dataset checks.

The official DLPFC donor names and retained run aliases are recorded in [dlpfc_sample_aliases.json](dlpfc_sample_aliases.json). These naming corrections do not change sample IDs, input bytes, sampling or saved results.
