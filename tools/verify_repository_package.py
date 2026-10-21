"""Check package bytes and local documentation without downloading or solving."""
from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote


ROOT = Path(__file__).resolve().parents[1]


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def inspect_package() -> dict:
    errors = []
    manifest = subprocess.run(
        [sys.executable, '-I', '-X', 'utf8', '-B', str(ROOT / 'tools/verify_manifest.py'), '--strict'],
        cwd=ROOT, capture_output=True, text=True, encoding='utf-8')
    if manifest.returncode:
        errors.append(manifest.stdout + manifest.stderr)

    binding = json.loads((ROOT / 'paper_protocol/source_binding.json').read_text(encoding='utf-8'))
    for name, wanted in binding['source_files'].items():
        if sha(ROOT / name) != wanted:
            errors.append('Bound source mismatch: ' + name)

    version = json.loads((ROOT / 'LOCAL_VERSION.json').read_text(encoding='utf-8'))
    packaged = version['standalone_repository']['bundled_file_sha256']
    for name, wanted in packaged.items():
        if not (ROOT / name).is_file() or sha(ROOT / name) != wanted:
            errors.append('Bundled attachment/resource mismatch: ' + name)

    input_manifest = json.loads((ROOT / 'data/input_files.json').read_text(encoding='utf-8'))
    for row in input_manifest['files']:
        if row.get('included_in_release'):
            path = ROOT / row['path']
            if not path.is_file() or path.stat().st_size != row['size_bytes'] or sha(path) != row['sha256']:
                errors.append('Bundled input-manifest mismatch: ' + row['path'])

    s2 = ROOT / 'lrot_output/s23_support/Supplementary_File_S2_synthetic_seed_pairs.csv'
    with s2.open(encoding='utf-8-sig', newline='') as handle:
        pairs = list(csv.DictReader(handle))
    if len(pairs) != 30 or [(r['seed_A'], r['seed_B']) for r in pairs] != [(str(1000+i), str(2000+i)) for i in range(30)]:
        errors.append('Selected S2 ordered seed-pair registry differs')

    links = []
    evidence_links = []
    for path in sorted(ROOT.rglob('*.md')):
        rel = path.relative_to(ROOT).as_posix()
        is_history = rel.startswith('provenance/') or '/history/' in '/' + rel
        for target in re.findall(r'!?\[[^\]]*\]\(([^\s)]+)(?:\s+[^)]*)?\)', path.read_text(encoding='utf-8-sig')):
            if re.match(r'^[a-z][a-z0-9+.-]*:', target, re.I) or target.startswith('#'):
                continue
            local = (path.parent / unquote(target.split('#')[0])).resolve()
            ok = local.is_relative_to(ROOT.resolve()) and local.exists()
            entry = {'document': rel, 'target': target, 'pass': ok}
            if is_history:
                evidence_links.append(entry)
            else:
                links.append(entry)
                if not ok:
                    errors.append('Current documentation link missing or outside repository: ' + rel + ' -> ' + target)

    review_hits = []
    for path in ROOT.rglob('*'):
        if path.is_file() and path.name.lower().startswith('readme') and path.suffix.lower() == '.md':
            if re.search(r'supervisor|mentor|导师|待审|审阅包|review package|review round', path.read_text(encoding='utf-8-sig'), re.I):
                review_hits.append(path.relative_to(ROOT).as_posix())
    if review_hits:
        errors.append('Internal-delivery wording in README: ' + ', '.join(review_hits))

    large_files = [{'path': p.relative_to(ROOT).as_posix(), 'size_bytes': p.stat().st_size}
                   for p in sorted(ROOT.rglob('*'))
                   if p.is_file() and '.git' not in p.relative_to(ROOT).parts and p.stat().st_size > 25 * 1024 * 1024]
    if not re.search(r'^\*\s+-text\s*$', (ROOT / '.gitattributes').read_text(encoding='utf-8'), re.M):
        errors.append('Exact-byte Git attribute rule is missing')
    return {
        'passed': not errors, 'errors': errors,
        'manifest_output': manifest.stdout.strip(), 'source_binding_files_checked': len(binding['source_files']),
        'bundled_resources_and_attachments_checked': len(packaged), 'S2_ordered_pairs_checked': len(pairs),
        'current_documentation_links_checked': len(links),
        'historical_documentation_link_exceptions': [r for r in evidence_links if not r['pass']],
        'files_above_browser_25_MiB': large_files,
        'raw_tissue_and_CellChatDB_inputs_bundled': False,
        'scope': 'Package bytes, local path references and selected metadata only. No network, array scoring, drawing, training, solving or whole-paper reproduction.'
    }


if __name__ == '__main__':
    report = inspect_package()
    print(json.dumps(report, ensure_ascii=False, indent=2))
    raise SystemExit(0 if report['passed'] else 1)
