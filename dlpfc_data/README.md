# DLPFC raw inputs

Supply expression H5 and tissue-position TXT files for 151507/151508/151509 (Br5292) and 151669/151670/151671 (Br5595), plus `barcode_level_layer_map.tsv`. File names, direct access patterns and checksum scope are specified in [data/README.md](../data/README.md) and `data/input_files.json`.

Verify before running:

```text
python tools/verify_inputs.py --dataset dlpfc --data-root "<parent-of-dlpfc_data>"
```

F5 accepts `--data-dir` pointing directly to this data directory. Other legacy loaders may use `LROT_DLPFC_DATA_DIR`; consult the selected entrypoint. The optional download helper does not verify historical identity and must be followed by checksum checks. Raw files are supplied separately, not deleted or replaced by packaging.

Historical records and unchanged code identifiers may use `Br8325` for the selected Br5292 sections and `Br6423` for the selected Br5595 sections. These are retained run aliases, not additional donors. See [the explicit mapping](../data/dlpfc_sample_aliases.json). Input hashes and section IDs, not a name replacement, bind the actual samples.
