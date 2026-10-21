"""Verbatim selected definitions from the hash-bound run_driver.py.
Configured only by paper_protocol.binding; no historical driver main is run.
"""
from pathlib import Path
import json,hashlib,inspect,time
import numpy as np

class RejectedPlan(RuntimeError):pass

def arrayhash(x):
    x=np.ascontiguousarray(x);h=hashlib.sha256();h.update(str(x.shape).encode());h.update(str(x.dtype).encode());h.update(x.tobytes());return h.hexdigest()

def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [clean(v) for v in x]
    if isinstance(x,np.generic):return clean(x.item())
    if isinstance(x,np.ndarray):return clean(x.tolist()) if x.size<=5000 else {'shape':list(x.shape),'dtype':str(x.dtype),'sha256':arrayhash(x)}
    if isinstance(x,float) and not np.isfinite(x):return None
    if x is None or isinstance(x,(str,int,float,bool)):return x
    return str(x)

def context():
    out={}
    for f in inspect.stack()[2:6]:
        if str(C) not in f.filename:continue
        v={k:clean(x) for k,x in f.frame.f_locals.items() if k in ['seed','ds','s','sigma','sig','mode','name','variant','gamma','g','alpha','beta','i','sample_A','sample_B','lr_program','eps_mult','method','rotation','n_spots'] and isinstance(x,(str,int,float,bool,np.generic))}
        out[Path(f.filename).name+':'+f.function]=v
    return out

def evaluate(P,A,D):
    Q=P/P.sum(1,keepdims=True);out={'entropy':float(-(Q*np.log(np.maximum(Q,1e-300))).sum(1).mean())}
    if 'region_labels' in A and 'region_labels' in D and max(np.max(A['region_labels']),np.max(D['region_labels']))>0:
        labs=np.asarray(D['region_labels'],int);la=np.asarray(A['region_labels'],int);n=int(max(labs.max(),la.max()))+1
        pred=(Q@np.eye(n)[labs]).argmax(1)
        from sklearn.metrics import adjusted_rand_score
        out.update(soft_acc=float((pred==la).mean()),soft_ari=float(adjusted_rand_score(la,pred)),same_mass=float(P[la[:,None]==labs[None,:]].sum()/P.sum()))
    else:out['annotation_metrics_available']=False
    return out

def solver(*a,**kw):
    bound=signature.bind(*a,**kw);bound.apply_defaults();v=bound.arguments
    A=v.pop('slice_A',None);D=v.pop('slice_B',None);S=v.pop('S_lr',None)
    # Resolve actual argument names from signature rather than assuming a name.
    if A is None:
        values=list(bound.arguments.values());A,D,S=values[:3];v={k:x for k,x in list(bound.arguments.items())[3:]}
    params=dict(v);params.update(strict=False,return_diagnostics=True,verbose=False,inner_max_iter=10000,marginal_abs_tol=1e-8,marginal_rel_tol=1e-5,tol=1e-5)
    identities={'A':{k:arrayhash(A[k]) for k in ['coords','expr','region_labels']},'B':{k:arrayhash(D[k]) for k in ['coords','expr','region_labels']},'gene_names_A':list(A['gene_names']),'gene_names_B':list(D['gene_names']),'lr':arrayhash(S),'params':clean(params),'solver':fingerprint}
    key=hashlib.sha256(json.dumps(identities,sort_keys=True).encode()).hexdigest()
    if driver=='serial_runtime_audit.py':key+='__trial_'+str(len(calls)+1)
    j=E/(key+'.json');p=E/(key+'.npz')
    ctx=context();calls.append({'key':key,'context':ctx,'call':len(calls)+1});(O/'calls.json').write_text(json.dumps(calls,ensure_ascii=False,indent=2),encoding='utf-8')
    if j.exists():
        record=json.loads(j.read_text(encoding='utf-8'))
        if record['accepted']:
            with np.load(p) as z:return z['P'],z['loss'].tolist()
        raise RejectedPlan('Previously unaccepted result; no silent replacement '+key)
    # Reuse only an exactly identical, already accepted configuration in this
    # isolated revalidation. Never reuse plans for timing/scaling measurements.
    if driver not in ['lrot_scaling.py','run_runtime_benchmark.py','serial_runtime_audit.py']:
        existing=next((x for x in (B/'runs').glob('*/evidence/'+key+'.json') if x!=j),None)
        if existing is not None:
            record=json.loads(existing.read_text(encoding='utf-8'))
            if record.get('accepted') and record['input_identity']==identities:
                import shutil
                record['exact_checkpoint_reuse']=str(existing);record['driver']=driver;record['context']=ctx
                j.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8');shutil.copy2(existing.with_suffix('.npz'),p)
                with np.load(p) as z:return z['P'],z['loss'].tolist()
    cert.reset_dual_cache();t=time.time();attempts=[]
    P,loss,d=actual(A,D,S,**params);attempts.append(d)
    if not d['success'] and d['stop_reason']=='outer_iteration_cap':
        params['max_iter']=1000;cert.reset_dual_cache();P,loss,d=actual(A,D,S,**params);attempts.append(d)
    core._last_compute_seconds=time.time()-t
    record={'key':key,'driver':driver,'stage':stage,'accepted':bool(d['success']),'input_identity':identities,'final_parameters':clean(params),'context':ctx,'runtime_seconds':core._last_compute_seconds,'diagnostics':d,'attempts':attempts,'metrics':evaluate(P,A,D)}
    j.write_text(json.dumps(clean(record),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8');np.savez_compressed(p,P=P,loss=np.array(loss))
    print('AUDIT',group,len(calls),d['success'],d['stop_reason'],round(record['runtime_seconds'],2),record['metrics'],flush=True)
    if not d['success']:raise RejectedPlan('Unaccepted plan cannot enter figure/report: '+key)
    return P,loss
