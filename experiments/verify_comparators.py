"""Independently check saved comparator plans, DLPFC embeddings and GO statistics."""
from pathlib import Path
import argparse,csv,hashlib,json,socket
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output-dir',required=True,type=Path)
a=p.parse_args();ROOT=Path(__file__).resolve().parents[1];OUT=a.output_dir.resolve()
if OUT.exists() or OUT==ROOT or ROOT in OUT.parents:p.error('Choose a new directory outside the package.')
OUT.mkdir(parents=True)
def deny(*a,**k):raise RuntimeError('Local verification; networking disabled')
socket.create_connection=deny;socket.socket.connect=deny;socket.socket.connect_ex=deny
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import hypergeom
from sklearn.metrics import adjusted_rand_score
from statsmodels.stats.multitest import multipletests
def arrayhash(x):
 x=np.ascontiguousarray(x);return hashlib.sha256(str(x.shape).encode()+str(x.dtype).encode()+x.tobytes()).hexdigest()
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
base=ROOT/'provenance/comparator_audit'
records=json.loads((base/'experiments.json').read_text('utf-8'))
plans=np.load(base/'plans.npz',allow_pickle=False);labs=np.load(base/'ordered_labels.npz',allow_pickle=False)
for r in records:
 P=plans[r['plan_key']];assert arrayhash(P)==r['plan_key'][2:]
 A=labs[r['dataset']+'_A'];B=labs[r['dataset']+'_B'];classes=np.eye(int(max(A.max(),B.max()))+1)[B]
 pred=(P@classes).argmax(1)
 assert float((pred==A).mean())==r['accuracy'] and float(adjusted_rand_score(A,pred))==r['ARI']
 assert r['input_hashes_captured_in_record'] or r['matching_input_bound_plan_available']
neural=[];base=ROOT/'provenance/dlpfc_comparators'
for r in json.loads((base/'results.json').read_text('utf-8')):
 d=base/r['pair'];la=np.load(d/'ordered_labels.npz',allow_pickle=False);train=np.load(d/'actual_training_input.npz',allow_pickle=False)
 assert sha(d/'actual_training_input.npz')==r['actual_training_input_sha256']
 assert sha(d/'embedding.npy')==r['embedding_sha256']
 assert train['X'].shape==(2000,50) and train['edgeList'].shape==(2,24000)
 rows=list(csv.DictReader((d/'ordered_spots.csv').open(encoding='utf-8')))
 assert len(rows)==2000 and all(x['batch']=='A' for x in rows[:1000]) and all(x['batch']=='B' for x in rows[1000:])
 assert np.array_equal(np.array([int(x['label']) for x in rows[:1000]]),la['A_labels'])
 assert np.array_equal(np.array([int(x['label']) for x in rows[1000:]]),la['B_labels'])
 logs=json.loads((d/'training_losses.json').read_text('utf-8'));assert len(logs)==1000 and sum(x['phase']=='joint' for x in logs)==500
 for kind,f in [('fresh','embedding.npy'),('historical','historical_embedding.npy')]:
  e=np.load(d/f,allow_pickle=False);nearest=cdist(e[:1000],e[1000:]).argmin(1)
  alternate=((e[:1000].astype(float)[:,None,:]-e[1000:].astype(float)[None,:,:])**2).sum(2).argmin(1)
  assert np.array_equal(nearest,alternate)
  pred=la['B_labels'][nearest];acc=float((pred==la['A_labels']).mean());ari=float(adjusted_rand_score(la['A_labels'],pred))
  ref=r if kind=='fresh' else json.loads((d/'historical_verification.json').read_text('utf-8'))
  assert acc==ref['accuracy'] and ari==ref['ARI']
  neural.append({'pair':r['pair'],'kind':kind,'accuracy':acc,'ARI':ari})
base=ROOT/'provenance/go_enrichment'
identity=json.loads((base/'resource_identity.json').read_text('utf-8'))
rows=list(csv.DictReader((base/'all_terms.csv').open(encoding='utf-8')))
assert len(rows)==14584
fg=np.array([int(r['a']) for r in rows]);bg=np.array([int(r['c']) for r in rows])
ps=hypergeom.sf(fg-1,identity['n_background'],bg,identity['n_foreground']);qs=multipletests(ps,method='fdr_bh')[1]
assert np.allclose(ps,[float(r['p']) for r in rows],rtol=1e-10,atol=0)
assert np.allclose(qs,[float(r['padj']) for r in rows],rtol=1e-10,atol=0)
display=json.loads((ROOT/'figures/supplement/inputs/S3_go_display_rows.json').read_text('utf-8'))['rows'];lookup={r['go']:r for r in rows}
for r in display:
 s=lookup[r['go']];assert r['term']==s['term'] and r['a']==int(s['a']) and r['c']==int(s['c']) and r['padj']==float(s['padj'])
result={'success':True,'comparator_plans_rescored':len(records),'DLPFC_scores':neural,'GO_terms_independently_checked':len(rows),'S3_display_rows_match':True,'training_executed':False,'whole_paper_certified':False}
(OUT/'verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n','utf-8');print(json.dumps(result))
