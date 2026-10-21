# Package contents

This directory is the repository root. `RELEASE_FILES.txt` is the authoritative
allowlist; `MANIFEST.sha256` and `MANIFEST_size.txt` verify the exact package.

Included: the original implementation and bound numerical protocol, selected
images and drawing inputs, retained results and provenance, fixed GO/MGI files
in [data/go](../data/go/README.md), and the unchanged Supplementary File S1
[workbook](../supplementary_files/LROT_Supplementary_Tables_S11_S12.xlsx).
Supplementary File S2 is [one CSV containing the selected 30 ordered seed pairs
and results](../lrot_output/s23_support/Supplementary_File_S2_synthetic_seed_pairs.csv).
The upstream MIT STAligner source recipe, licence and disclosed patch remain
included. Resource-specific attribution and terms remain separate from project MIT.

Not included: original tissue matrices/coordinates/annotations, CellChatDB
binaries, dependency environments, restricted fonts, manuscript DOCX/PDF files,
credentials, private backups or the full local experiment collection.
Use [input preparation](../data/README.md) and [external resources](external_resources.md).
Missing external inputs must fail verification; no substitute data are supplied.

Run `python tools/verify_repository_package.py` for package/preparation checks.
Build another exact directory using `python tools/build_release.py --output
"<new-directory>"`, then run its `python tools/verify_manifest.py --strict`.
The output must be outside the source tree and must not already exist.
See [upload/layout instructions](repository_upload.md).

Package integrity, saved-image export, saved-result drawing, independent scoring
and end-to-end experiment reproduction remain distinct. Case-specific evidence
and unresolved scientific limits are preserved in [the verification scope](release_readiness.md).
Historical records may retain earlier names, file counts and directory layouts;
they are not the current package inventory. No new software version or DOI is assigned.
