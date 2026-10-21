"""New explicitly defined full STAligner baseline, not a historical-driver reproduction.
Spatial six-nearest-neighbour edges within each section; PCA-50 features;
package v1.0.0, CPU, Adam .001, seed666, n_epochs1000 = 500 pretrain +500 joint.
The archived n_epochs200 call never entered joint alignment and is excluded.
"""
from pathlib import Path
import os,json,time,numpy as np
from scipy import sparse
from scipy.spatial.distance import cdist
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import adjusted_rand_score
from real_st_loader import RealisticSTGenerator
from run_dlpfc_baselines import load_dlpfc_slice,apply_misalignment
import scanpy as sc
import STAligner as sta
OUT=Path(os.environ['LROT_OUTPUT_DIR']);OUT.mkdir(exist_ok=True,parents=True)
def run_pair(name,A,B):
    start=time.time();ads=[];edges=[];offset=0
    for batch,sd in enumerate([A,B]):
        ad=sc.AnnData(X=sd['expr'].astype(np.float32));ad.var_names=list(sd['gene_names']);ad.obs_names=[f'{batch}_{i}' for i in range(sd['n_spots'])]
        ad.obs['batch_name']=str(batch);ad.obsm['spatial']=sd['coords'].astype(np.float32)
        nbr=NearestNeighbors(n_neighbors=7).fit(sd['coords']).kneighbors(sd['coords'],return_distance=False)[:,1:]
        u=np.repeat(np.arange(sd['n_spots']),6);v=nbr.ravel()
        edge=np.concatenate([np.array([u,v]),np.array([v,u])],axis=1)+offset
        edges.append(edge);offset+=sd['n_spots'];ads.append(ad)
    import anndata as an
    ad=an.concat(ads,join='inner',merge='same');ad.obs['batch_name']=ad.obs['batch_name'].astype('category')
    sc.pp.normalize_total(ad,target_sum=1e4);sc.pp.log1p(ad);sc.tl.pca(ad,n_comps=50)
    train=sc.AnnData(X=sparse.csr_matrix(ad.obsm['X_pca'].astype(np.float32)));train.obs=ad.obs.copy();train.uns['edgeList']=np.concatenate(edges,axis=1)
    sta.train_STAligner(train,n_epochs=1000,lr=.001,verbose=False,device='cpu',random_seed=666,knn_neigh=50)
    X=train.obsm['STAligner'];n=A['n_spots'];D=cdist(X[:n],X[n:]);pred=B['region_labels'][D.argmin(1)]
    assert np.isfinite(X).all()
    row={'dataset':name,'accuracy_nn':float(np.mean(pred==A['region_labels'])),'ari_nn':float(adjusted_rand_score(A['region_labels'],pred)),'end_to_end_seconds':time.time()-start,'spot_counts':[A['n_spots'],B['n_spots']],'graph_edges':int(train.uns['edgeList'].shape[1]),'protocol':'PCA50 spatial KNN6 undirected within-section graphs; Adam .001; seed666; 500 pretraining +500 joint epochs; MNN50','interpretation':'new complete baseline; excludes the pretraining-only archived adaptation; latent NN metrics are not transport entropy or calibration'}
    print('FULL_STALIGNER',row,flush=True);return row
def main():
    rows=[]
    gen=RealisticSTGenerator(n_spots_A=1000,n_spots_B=1000,n_genes=150,seed=42);p=gen.generate_paired_slices(rotation=5.,batch_effect=.15)
    rows.append(run_pair('realistic-seed42',p.slice_A.to_dict(),p.slice_B.to_dict()));(OUT/'full_protocol_results.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
    for a,b in [('151507','151508'),('151508','151509')]:
        A=load_dlpfc_slice(a,1000,42);B=load_dlpfc_slice(b,1000,99);B['coords']=apply_misalignment(B['coords'],rotation_deg=2.,translation=.5)
        rows.append(run_pair(a+'-'+b,A,B));(OUT/'full_protocol_results.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
if __name__=='__main__':main()
