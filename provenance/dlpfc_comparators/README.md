# DLPFC comparator verification

The saved historical embeddings independently rescore to accuracy .580/.518 and ARI .383938/.315308 in reconstructed section/barcode order. Historical incoming training-input hashes were not captured; this fact is retained.

Current Tables S9/S22 use new explicitly input-bound, single-thread deterministic runs: accuracy .586/.525 and ARI .377807/.302971. Each has 500 STAGATE pretraining and 500 STAligner joint steps, seed 666, PCA50, within-section KNN6 plus reverse edges and MNN50. Each pair uses source/target sampling seeds 42/99; the shared section is not necessarily the same subset in both roles. Labels are withheld from training and used only for independent A-to-B latent-nearest-neighbour scoring. These are two prescribed technical runs, not independent donors or a stability distribution. Their ranking against the reported comparator accuracies is unchanged.

`actual_training_input.npz`, `ordered_spots.csv`, raw-file/array identities, loss logs, final model weights and returned embeddings bind the new computations. Preprocessing and training are separate processes to release dense transcriptome memory. The package-returned embedding is retained without redefining it as an extra post-update evaluation.

```text
python experiments/prepare_dlpfc_staligner.py --data-root "<data-root>" --output-dir "<new-prepared>"
python tools/prepare_external_resources.py staligner --output-dir "<new-fixed-source>"
python experiments/train_dlpfc_staligner.py --input-dir "<new-prepared>" --source-dir "<new-fixed-source>" --output-dir "<new-trained>"
python experiments/verify_comparators.py --output-dir "<new-check>"
```

The same pinned scientific dependency versions as F3 are required. This verification does not certify every manuscript experiment or cross-platform repeatability.
