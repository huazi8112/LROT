"""Read fixed sampled DLPFC rows while retaining barcode and layer order."""
from pathlib import Path
import hashlib
import numpy as np
import scanpy as sc
layers_order=['L1','L2','L3','L4','L5','L6','WM']
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def arrayhash(x):
 x=np.ascontiguousarray(x);h=hashlib.sha256(str(x.shape).encode()+str(x.dtype).encode());h.update(x.tobytes());return h.hexdigest()
def load_dlpfc_input(data_root, sample, seed):
    data = Path(data_root)/'dlpfc_data'
    paths = [data/f'{sample}_filtered_feature_bc_matrix.h5',
             data/f'{sample}_tissue_positions_list.txt', data/'barcode_level_layer_map.tsv']
    ad = sc.read_10x_h5(paths[0])
    if not ad.var_names.is_unique:
        ad.var_names_make_unique()
    pos = {}
    for line in paths[1].read_text(encoding='utf-8').splitlines():
        c = line.split(',')
        if len(c) == 6:
            pos[c[0]] = (float(c[4]), float(c[5]))
    labels = {}
    for line in paths[2].read_text(encoding='utf-8').splitlines():
        c = line.split('\t')
        if len(c) == 3 and c[1] == sample:
            labels[c[0]] = c[2]
    barcodes = [b for b in ad.obs_names if b in pos and b in labels]
    coords = np.array([pos[b] for b in barcodes], dtype=np.float64)
    coords -= coords.mean(axis=0)
    idx = np.random.RandomState(seed).choice(len(barcodes), min(1000, len(barcodes)), replace=False)
    selected = [barcodes[i] for i in idx]
    sd = dict(expr=ad[selected].X.toarray().astype(np.float32), coords=coords[idx],
              region_labels=np.array([layers_order.index(labels[b]) for b in selected], dtype=int),
              gene_names=list(ad.var_names), n_spots=len(idx), sample_id=sample,
              barcodes=selected, sampled_indices=idx)
    identity = {'sample':sample, 'seed':seed, 'spots':len(idx), 'eligible_spots':len(barcodes),
                'raw_files':{p.name:sha(p) for p in paths},
                'arrays':{k:arrayhash(sd[k]) for k in ['expr','coords','region_labels']},
                'gene_names_sha256':arrayhash(np.array(sd['gene_names'],dtype=str)),
                'barcodes_sha256':arrayhash(np.array(selected,dtype=str))}
    return sd, identity
