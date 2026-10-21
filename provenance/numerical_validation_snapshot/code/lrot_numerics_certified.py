"""Optional difficult-temperature solver using the SAME balanced OT problem.

Residual acceptance, never optimizer status alone, determines success.
Evaluated separately from the primary numerical revision; archived failures
remain available. No modification to LR/fusion costs or regularization.
"""
import numpy as np
from scipy.special import logsumexp
from scipy.optimize import minimize
from lrot_numerics import balanced_sinkhorn as scaling, plan_diagnostics, TransportConvergenceError

_DUAL_CACHE = {}
def reset_dual_cache():
    _DUAL_CACHE.clear()

def feasible_round(P,a,b):
    R = P * np.minimum(a/P.sum(1),1)[:,None]
    R = R * np.minimum(b/R.sum(0),1)[None,:]
    g = np.maximum(a-R.sum(1),0);h = np.maximum(b-R.sum(0),0)
    if g.sum() > 0:
        R = R + np.outer(g,h)/g.sum()
    return R

def balanced_sinkhorn(log_kernel,a,b,max_iter=10000,abs_tol=1e-8,rel_tol=1e-5,strict=True,check_every=10):
    P,info=scaling(log_kernel,a,b,max_iter=min(max_iter,1000),abs_tol=abs_tol,rel_tol=rel_tol,strict=False,check_every=check_every)
    if info['converged']:
        info['strategy']='scaling_then_dual_if_needed'
        return P,info
    L=np.asarray(log_kernel,dtype=float);a=np.asarray(a);b=np.asarray(b)
    L=L-L.max(1,keepdims=True);L=L-L.max(0,keepdims=True)
    # Warm-start dual potentials with unprotected exact scaling when representable.
    if L.min()>-650:
        K=np.exp(L);u=a.copy();v=b.copy()
        for it in range(500):
            v=b/(K.T@u);u=a/(K@v)
        initial=np.log(v)
    else:
        lu=np.log(a);initial=np.log(b)
        for it in range(300):
            initial=np.log(b)-logsumexp(L+lu[:,None],axis=0)
            lu=np.log(a)-logsumexp(L+initial[None,:],axis=1)
    if len(b) in _DUAL_CACHE:
        initial = _DUAL_CACHE[len(b)].copy()
    initial=initial-initial[-1]
    scale=1/float(b.min());loga=np.log(a)
    def fg(x):
        v=np.r_[x,0.];Z=L+v[None,:];lse=logsumexp(Z,axis=1)
        p=np.exp(loga[:,None]+Z-lse[:,None])
        f=scale*(a@lse-b@v)
        g=scale*(p.sum(0)-b)
        return float(f),g[:-1]
    opt=minimize(fg,initial[:-1],method='L-BFGS-B',jac=True,
        options={'gtol':1e-10,'ftol':1e-16,'maxiter':3000,'maxls':50,'maxcor':30})
    v=np.r_[opt.x,0.];Z=L+v[None,:];lse=logsumexp(Z,axis=1)
    P=np.exp(loga[:,None]+Z-lse[:,None])
    _DUAL_CACHE[len(b)] = v.copy()
    raw = P.copy()
    rounded = feasible_round(raw,a,b)
    correction_l1 = float(abs(rounded-raw).sum())
    positive = rounded > 0
    primal = float(-(L*rounded).sum() + (rounded[positive]*(np.log(rounded[positive])-1)).sum())
    dual = float(a@(loga-lse)+b@v-raw.sum())
    gap = primal-dual
    # A feasibility correction is allowed only with a measured small mass
    # perturbation AND a convex-dual upper bound on subproblem suboptimality.
    certificate = bool(correction_l1 <= 1e-6 and -1e-8 <= gap <= 1e-5)
    raw_diag = plan_diagnostics(raw,a,b)
    if certificate:
        P = rounded
    diag=plan_diagnostics(P,a,b)
    ok=diag['finite'] and diag['nonnegative'] and max(diag['row_max_abs'],diag['col_max_abs'])<=abs_tol and max(diag['row_max_relative'],diag['col_max_relative'])<=rel_tol
    out=dict(iterations=info['iterations']+int(opt.nit),converged=bool(ok),stop_reason='marginal_tolerance' if ok else 'dual_polish_residual_failure',algorithm='reduced_balanced_dual_LBFGS',strategy='scaling_dual_and_bounded_feasibility_correction',optimizer_success=bool(opt.success),optimizer_message=str(opt.message),optimizer_iterations=int(opt.nit),initial_scaling_iterations=info['iterations'],raw_marginal_diagnostics=raw_diag,feasibility_correction_used=certificate,feasibility_correction_l1=correction_l1,dimensionless_primal_dual_gap=gap,duality_gap_limit=1e-5,correction_l1_limit=1e-6,abs_tol=abs_tol,rel_tol=rel_tol,**diag)
    if not ok and strict:raise TransportConvergenceError('Dual polishing failed numerical acceptance',out)
    return P,out
