# Tested clean-environment installation

A new Python 3.10.11 venv was created on Windows AMD64 with `include-system-site-packages = false`. User site packages were disabled. Dependencies were downloaded from the authorized public PyPI/PyG sources; none were copied from the original research installation. Scientific module origins were checked inside the new venv. The original research environment was not modified.

## Minimal plotting installation

From the repository root in PowerShell, use Python 3.10.11 and new external directories:

```powershell
python --version
python -m venv ../lrot_plot_env
../lrot_plot_env/Scripts/python.exe -m pip install --index-url https://pypi.org/simple --only-binary=:all: -r docs/environments/plotting.txt
../lrot_plot_env/Scripts/python.exe tools/verify_manifest.py
../lrot_plot_env/Scripts/python.exe figures/render_current.py --figure F2 --mode records --output-dir ../lrot_demo_f2
```

The clean-install F2 PNG matched the approved asset byte for byte. The example uses retained inputs and does not execute a solver. The original README's equivalent three-package command was also tested. The installation is performed by these explicit pip commands, not by a numerical entrypoint.

## Fixed F3 training installation

Create a separate environment or extend the plotting environment with the [fixed Windows CPU profile](environments/f3-windows.txt). This profile includes the plotting packages, core numerical packages, fixed training dependencies and `rdata` for the resource preparation tool.

```powershell
python -m venv ../lrot_f3_env
../lrot_f3_env/Scripts/python.exe -m pip install --index-url https://pypi.org/simple --only-binary=:all: -r docs/environments/f3-windows.txt
../lrot_f3_env/Scripts/python.exe -m pip check
../lrot_f3_env/Scripts/python.exe tools/prepare_external_resources.py staligner --output-dir ../lrot_external/staligner
$env:PYTHONPATH = (Resolve-Path ../lrot_external/staligner).Path
../lrot_f3_env/Scripts/python.exe experiments/reproduce_f3_staligner.py --mode train --output-dir ../lrot_f3_clean_check
```

The profile supplies `torch==2.11.0` from the Windows PyPI wheel and the matching `torch-scatter==2.1.2+pt211cpu` / `torch-sparse==0.6.18+pt211cpu` CPU wheels from the fixed PyG wheel page. These two compiled extensions were missing from the earlier dependency inventory and are now listed explicitly. STAligner itself is reconstructed offline from the bundled fixed upstream sources and patch; it is not fetched as an arbitrary PyPI package. Keep annoy/hnswlib absent for this fixed exact-neighbour protocol. All original training gates remain active.

The [complete tested dependency lock](environments/verified-windows-lock.txt) records the resolved transitive versions. It can replace the F3 profile in the installation command when the same Python/Windows CPU configuration is required. The bootstrap pip version is an installation tool, not included as a scientific dependency; the separately prepared STAligner source is also excluded from the pip lock. The broad root `requirements.txt` remains a research inventory; the clean-install result does not certify all optional baseline/document-processing packages in that inventory.

If a proxy is required, add `--proxy "<proxy-url>"` to the pip command using your own network configuration. The verification used a process-level proxy, a project-local cache and isolated pip configuration. TLS verification remained enabled; no system settings were changed.

## Results and boundaries

- F2 saved-result drawing: PNG byte-identical to the approved figure.
- F2 fresh diagnostics: all four 80-update histories matched within 2.22×10⁻¹⁶; the retained 123-value main trace was unchanged and is not claimed newly executed.
- F3 saved-embedding scoring: both retained embeddings produced agreement .930 and ARI .8535028559314339.
- F3 new full training: 1000 steps in the clean dependency environment, agreement .930 and the same ARI. Actual training input, returned embedding, all losses and all final weight tensors matched the recorded fixed case. Seed and protocol were unchanged.
- Dependency checks passed. The full training entrypoint blocks network connections.

Evidence: [verification.json](../provenance/clean_environment/verification.json) and [dependency versions/download hashes](../provenance/clean_environment/installed_dependencies.json). Detailed installation reports, full console output and the environment itself remain in the private local audit outside the code distribution.

This validates the listed tasks on the same Windows computer. It does not establish a new-machine/cross-platform result, recover historical .934, certify every experiment, or authenticate earlier resource acquisition history. Manuscript/PDF finalization and author-approved publication metadata remain separate.
