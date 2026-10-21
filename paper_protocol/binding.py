"""Bind archived numerical functions without executing the historical run driver.
No original data search, locking loop, experiment main or temporary source patch.
The preserved wrapper uses only a NEW isolated output tree for any checkpoints.
"""
from pathlib import Path
import ast,hashlib,json,sys,inspect,importlib,os
ROOT=Path(__file__).resolve().parents[1]
SNAP=ROOT/'provenance/numerical_validation_snapshot'
def verify():
    m=json.loads((Path(__file__).with_name('source_binding.json')).read_text(encoding='utf8'))
    for rel,h in m['source_files'].items():
        if hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=h:raise RuntimeError('Source hash mismatch: '+rel)
    sources=[(SNAP/'portable_run_driver.py').read_text(encoding='utf8'),Path(__file__).with_name('archived_wrapper.py').read_text(encoding='utf8')]
    for src in sources:
        seg={n.name:ast.get_source_segment(src,n) for n in ast.parse(src).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
        for name,h in m['extracted_definitions'].items():
            if hashlib.sha256(seg[name].encode()).hexdigest()!=h:raise RuntimeError('Wrapper definition changed: '+name)
    return m
def bind(output_root):
    m=verify();out=Path(output_root).resolve()
    if out.exists() or out==ROOT or ROOT in out.parents:raise ValueError('A fresh output directory outside the bundle is required')
    os.environ.update(MPLBACKEND='Agg',MPLCONFIGDIR=str(out/'mpl_cache'),OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2')
    code=(SNAP/'code').resolve()
    for name in ['lrot_core','lrot_numerics','lrot_numerics_certified','real_st_loader']:
        if name in sys.modules and Path(sys.modules[name].__file__).resolve()!=code/(name+'.py'):
            raise RuntimeError('Conflicting module already imported: '+name)
    sys.path.insert(0,str(code))
    core=importlib.import_module('lrot_core');cert=importlib.import_module('lrot_numerics_certified')
    for name in ['lrot_core','lrot_numerics','lrot_numerics_certified']:
        if Path(sys.modules[name].__file__).resolve()!=code/(name+'.py'):raise RuntimeError('Wrong module origin')
    core.balanced_sinkhorn=cert.balanced_sinkhorn
    from . import archived_wrapper as w
    w.C=code;w.B=out;w.driver='paper_protocol.py';w.stage='';w.group='paper_protocol'
    w.O=out/'runs'/w.group;w.E=w.O/'evidence';w.E.mkdir(parents=True,exist_ok=False)
    w.calls=[];w.hashes={};w.core=core;w.cert=cert;w.actual=core.fgw_lr_solver
    from threadpoolctl import threadpool_limits
    w._thread_limiter=threadpool_limits(2)
    w.signature=inspect.signature(w.actual);w.fingerprint=m['solver_fingerprint']
    (out/'binding.json').write_text(json.dumps(m,indent=2),encoding='utf8')
    return w

def solve_locked(wrapper,A,B,S,expected_identity,**parameters):
    """Refuse input/protocol differences before calling the unchanged wrapper."""
    bound=wrapper.signature.bind(A,B,S,**parameters);bound.apply_defaults()
    vals=list(bound.arguments.items());params=dict(vals[3:])
    params.update(strict=False,return_diagnostics=True,verbose=False,inner_max_iter=10000,
                  marginal_abs_tol=1e-8,marginal_rel_tol=1e-5,tol=1e-5)
    actual={'A':{k:wrapper.arrayhash(A[k]) for k in ['coords','expr','region_labels']},
            'B':{k:wrapper.arrayhash(B[k]) for k in ['coords','expr','region_labels']},
            'gene_names_A':list(A['gene_names']),'gene_names_B':list(B['gene_names']),
            'lr':wrapper.arrayhash(S),'params':wrapper.clean(params),'solver':wrapper.fingerprint}
    if actual!=expected_identity:raise ValueError('Input or protocol identity mismatch; solver not called')
    (wrapper.O/'locked_input_identity.json').write_text(json.dumps(actual,ensure_ascii=False,indent=2),encoding='utf8')
    return wrapper.solver(A,B,S,**parameters)
