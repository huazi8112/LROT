# F2 four-curve reconstruction from current code

Status: completed as a new current-code calculation; not authenticated as the missing historical run.

- Source SHA-256: 4ea83368f0651476389e0a86874c236afabd098c62a1a0c700aa734f5da966b7
- Synthetic design SHA-256: 8ac624e43b55c60287ff484cef9afc0c2c28640bdd39adc666a01073817d3db8
- Existing archive SHA-256: ef1d468bc2c233b5916f71feb57b5443128f410f99f1fde0ecf81de4e8f5f323 (main loss history: 123 points; no four ablation histories)
- Input identity: coordinates, labels, feature positions, and expression values matched the archive at stored precision; the two archived symbols Artemin/Neurturin map to the current names Artn/Nrtn using the manifest's explicit aliases.
- Current code calls: gamma=0/0.05/0.1/0.2; alpha=.5, beta=.3, max_iter=80, tol=1e-5; Sinkhorn max_iter=200, tol=1e-4.
- Output is a new standalone F panel. Original approved image, manuscript, and experiment archive were not changed.
- It does not prove the historical plotted curves came from this exact source/environment; original per-run logs/curve arrays and historical source/environment hashes remain absent.

| gamma | updates | first loss | final loss | curve SHA-256 |
|---:|---:|---:|---:|---|
| 0.00 | 80 | 0.545996514242 | 0.515895842511 | edc151ac28a33df3fc0e83b912674eed432bfe179fbe4deab313cf85d6ea2dc2 |
| 0.05 | 80 | 0.556951054676 | 0.527130223424 | 7c3ec9a510cfe4aee545c8dae372b7ef97e56bb2662478b0ff6e106242371e5b |
| 0.10 | 80 | 0.567792622937 | 0.537738936734 | 15315b2928e9ea473f721f823395d92e5d85e581ebea61794a53311639d2f7b7 |
| 0.20 | 80 | 0.589139221053 | 0.558709282765 | e97fdc5a2f10e092da5208c01d5489ac2e52ea2d20f3653ed6243f2df68e1d7b |
