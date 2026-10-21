"""Prepare one prespecified prospective F7 input from the four pinned local raw files.

This script never invokes an OT solver and never writes historical result paths.
"""
from __future__ import annotations
import argparse, hashlib, json, os, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW_DEFAULT = ROOT / 'cancer_data'
REFERENCE = ROOT / 'provenance/f7_prospective/reference_input_identity.json'
HIST_SHA = '5599a3ca51f245ed16d646fd8b4a46170121091727f91c4a63cf1435c0633690'
SAMPLES = ('V1_Breast_Cancer_Block_A_Section_1', 'V1_Breast_Cancer_Block_A_Section_2')

def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def pinned_raw_hashes() -> dict[str, str]:
    result = {}
    for line in (ROOT/'cancer_data/BREAST_SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
        if line and not line.startswith('#'):
            digest, rel = line.split(None, 1)
            result[rel.strip().replace('\\', '/')] = digest.lower()
    if len(result) != 4:
        raise RuntimeError('Expected exactly four pinned breast raw files')
    return result

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--data-dir', type=Path, default=RAW_DEFAULT,
                   help='Directory containing the four raw files; defaults to repository cancer_data/')
    p.add_argument('--historical-npz', type=Path,
                   help='Optional local old full NPZ for an extra archive input comparison')
    args = p.parse_args()
    raw = args.data_dir.resolve()
    historical = args.historical_npz.resolve() if args.historical_npz else None
    reference = json.loads(REFERENCE.read_text(encoding='utf-8'))
    out = args.output_dir.resolve()
    if not out.is_dir() or any(out.iterdir()):
        raise RuntimeError('Output directory must exist and be empty')
    raw_hashes = pinned_raw_hashes()
    if raw_hashes != reference['raw_file_sha256']:
        raise RuntimeError('Repository raw hash list differs from frozen reference')
    for rel, expected in raw_hashes.items():
        path = raw / rel
        if not path.is_file() or sha(path) != expected:
            raise RuntimeError('Raw input absent or hash mismatch: '+str(path))
    if historical and (not historical.is_file() or sha(historical) != HIST_SHA):
        raise RuntimeError('Optional archived NPZ absent or hash mismatch')
    os.environ['LROT_DATA_DIR'] = str(raw)
    sys.path.insert(0, str(ROOT/'experiments'))
    import run_cancer_visium_truecoord as source
    from paper_protocol.archived_wrapper import arrayhash
    if (source.N_SAMPLE, source.SEED_A, source.SEED_B, source.GAMMA) != (3000, 42, 99, 0.1):
        raise RuntimeError('Source settings differ from frozen prospective plan')
    A = source.load_breast_slice(SAMPLES[0], 3000, 42)
    B = source.load_breast_slice(SAMPLES[1], 3000, 99)
    B['coords'] = source.apply_misalignment(B['coords'], rotation_deg=2.0, translation=0.5)
    matches = None
    if historical:
        matches = {}
        with np.load(historical, allow_pickle=True) as z:
            for name, section in [('A', A), ('B', B)]:
                matches[name] = {
                    'barcodes': np.array_equal(np.asarray(section['barcodes'], dtype=str), z[name+'_barcodes'].astype(str)),
                    'coordinates': np.array_equal(section['coords'], z[name+'_coords']),
                    'expression_float16_archive': np.array_equal(section['expr'].astype(np.float16), z[name+'_expr']),
                    'gene_names': np.array_equal(np.asarray(section['gene_names'], dtype=str), z[name+'_gene_names'].astype(str)),
                }
            if z['P'].shape != (3000, 3000) or z['P_fgw'].shape != (3000, 3000):
                raise RuntimeError('Archived plan shapes differ from fixed input')
        print('archive_input_match', json.dumps(matches, sort_keys=True), flush=True)
        if not all(all(d.values()) for d in matches.values()):
            raise RuntimeError('Raw-loader input is not identical to archived input; no solve is allowed')
    S, top = source.compute_lr_strength_human(A, B, source.HUMAN_LR_DB)
    if S.shape != (3000, 3000) or not np.isfinite(S).all():
        raise RuntimeError('LR matrix failed fixed shape/finite checks')
    prior_hash = sha(ROOT/'experiments/run_cancer_visium.py')
    inputs = {
        'A': {k: arrayhash(A[k]) for k in ['coords','expr','region_labels']},
        'B': {k: arrayhash(B[k]) for k in ['coords','expr','region_labels']},
        'S_lr': arrayhash(S),
    }
    if (sha(ROOT/'experiments/run_cancer_visium_truecoord.py') != reference['source_script_sha256']
            or prior_hash != reference['human_prior_code_sha256']
            or inputs != reference['ordered_input_arrayhash']):
        raise RuntimeError('Ordered raw input, LR matrix or source code differs from frozen reference')
    barcode_hashes = {name: hashlib.sha256('\n'.join(sec['barcodes']).encode()).hexdigest()
                      for name, sec in [('A', A), ('B', B)]}
    if barcode_hashes != reference['ordered_barcode_sha256_newline_utf8']:
        raise RuntimeError('Ordered barcode identity differs from frozen reference')
    prep = out/'prepared_f7_raw_input.npz'
    np.savez_compressed(prep,
        A_coords=A['coords'], A_expr=A['expr'], A_region_labels=A['region_labels'],
        A_gene_names=np.asarray(A['gene_names'], dtype=str), A_barcodes=np.asarray(A['barcodes'], dtype=str),
        B_coords=B['coords'], B_expr=B['expr'], B_region_labels=B['region_labels'],
        B_gene_names=np.asarray(B['gene_names'], dtype=str), B_barcodes=np.asarray(B['barcodes'], dtype=str),
        S_lr=S,
    )
    record = {
        'scope': 'prospective raw-input preparation only; no solver was called',
        'raw_file_sha256': raw_hashes,
        'historical_npz_path': str(historical) if historical else None,
        'historical_npz_sha256': HIST_SHA if historical else None,
        'reference_input_identity_sha256': sha(REFERENCE),
        'ordered_barcode_sha256_newline_utf8': barcode_hashes,
        'archived_input_comparison': matches,
        'ordered_input_arrayhash': inputs,
        'source_script_sha256': sha(ROOT/'experiments/run_cancer_visium_truecoord.py'),
        'human_prior_code_sha256': prior_hash,
        'preparation_script_sha256': sha(Path(__file__)),
        'sample_ids': list(SAMPLES), 'seeds': {'A':42,'B':99},
        'n_spots':3000, 'B_transform': {'rotation_degrees':2.0,'translation':0.5},
        'human_lr_pairs':len(source.HUMAN_LR_DB),
        'matched_lr_pairs':sum(1 for lig,rec in source.HUMAN_LR_DB if lig.upper() in {str(g).upper() for g in A['gene_names']} and rec.upper() in {str(g).upper() for g in B['gene_names']}),
        'prepared_npz_path':str(prep), 'prepared_npz_sha256':sha(prep),
        'top_pairs_preview':[[str(l),str(r),float(w)] for l,r,w in top[:10]],
    }
    (out/'input_identity.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('prepared_npz_sha256', record['prepared_npz_sha256'], flush=True)
    print('input_identity', str(out/'input_identity.json'), flush=True)

if __name__ == '__main__':
    main()
