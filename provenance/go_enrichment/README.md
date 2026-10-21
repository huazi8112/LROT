# GO resource and result verification

The fixed local OBO and MGI GAF independently regenerate all 14,584 term records, including the 15 Figure S3 display rows, exactly. Foreground/background sizes remain 109/27,431; 1,815 terms pass adjusted p < .05 and at least three foreground hits. Independent hypergeometric and statsmodels BH checks agree.

Official current downloads on 2026-10-16 are byte-identical to the fixed inputs. GO data/products and MGI data/annotations declare CC BY 4.0 on the linked official policies in `resource_identity.json`. Original acquisition dates are unknown; this verification date is not a download-history reconstruction. The OBO header `releases/2026-07-26` is retained as a header identifier: the guessed release-directory URLs returned 404, so they are not offered as valid access URLs.

Raw files remain outside the code package. A separate local `LROT_GO_resources.zip` companion freezes these exact public bytes with attribution and checksums for the author to archive with the release. Until that archive is published, current URLs may change; checksum mismatch must fail instead of substituting a newer resource.

```text
python experiments/go_enrichment.py --obo "<go-basic.obo>" --gaf "<mgi.gaf.gz>" --output-dir "<new-out>"
```

The historical `lrot_output/lrot_go_enrichment_results.txt` is preserved and is not the current paper analysis. Current results are `all_terms.csv` and `lrot_go_enrichment_results.txt` in this directory. Neither GO enrichment nor resource membership validates communication or alignment correctness.
