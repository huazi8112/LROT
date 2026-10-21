"""Local numerical revision. No network or file-writing side effects.

Balanced OT uses the same kernel exp(-cost/epsilon), with row/column gauge
shifts (which do not change the constrained optimum). Ordinary scaling has
no additive denominator. Underflow/overflow cases use log-domain updates.
"""
import numpy as np
from scipy.special import logsumexp

class TransportConvergenceError(RuntimeError):
    def __init__(self, message, diagnostics):
        super().__init__(message)
        self.diagnostics = diagnostics

def plan_diagnostics(P, a=None, b=None):
    P = np.asarray(P, dtype=float)
    if P.ndim != 2 or min(P.shape) == 0:
        raise ValueError('Transport plan must be a nonempty matrix')
    a = np.full(P.shape[0], 1 / P.shape[0]) if a is None else np.asarray(a)
    b = np.full(P.shape[1], 1 / P.shape[1]) if b is None else np.asarray(b)
    finite = bool(np.isfinite(P).all())
    if not finite:
        return dict(finite=False, nonnegative=False, mass=None, row_max_abs=None,
                    col_max_abs=None, row_max_relative=None, col_max_relative=None)
    r, c = P.sum(1), P.sum(0)
    return dict(finite=True, nonnegative=bool((P >= 0).all()), mass=float(P.sum()),
                row_max_abs=float(abs(r-a).max()), col_max_abs=float(abs(c-b).max()),
                row_max_relative=float((abs(r-a)/a).max()),
                col_max_relative=float((abs(c-b)/b).max()))

def validate_plan(P):
    P = np.asarray(P, dtype=float)
    if P.ndim != 2 or not np.isfinite(P).all() or (P < 0).any():
        raise ValueError('Plan must be finite, nonnegative and two-dimensional')
    if P.sum() <= 0 or (P.sum(1) <= 0).any():
        raise ValueError('Zero-mass transport plan/row is undefined for evaluation')
    return P

def balanced_sinkhorn(log_kernel, a, b, max_iter=10000, abs_tol=1e-8,
                      rel_tol=1e-5, strict=True, check_every=10):
    L = np.asarray(log_kernel, dtype=np.float64)
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if L.shape != (len(a), len(b)) or not np.isfinite(L).all():
        raise ValueError('A finite log kernel matching both marginal sizes is required')
    if not np.isfinite(a).all() or not np.isfinite(b).all() or (a<=0).any() or (b<=0).any():
        raise ValueError('Marginals must be finite and strictly positive')
    if not np.isclose(a.sum(), b.sum(), rtol=0, atol=1e-12):
        raise ValueError('Marginals must have equal total mass')
    if max_iter < 1 or abs_tol <= 0 or rel_tol <= 0 or check_every < 1:
        raise ValueError('Invalid iteration limit or tolerance')
    # These gauge shifts change the cost only by row/column constants.
    L = L - L.max(axis=1, keepdims=True)
    L = L - L.max(axis=0, keepdims=True)
    log_mode = L.min() < -650
    loga, logb = np.log(a), np.log(b)
    if log_mode:
        lu, lv = loga.copy(), logb.copy()
    else:
        K = np.exp(L)
        u, v = a.copy(), b.copy()
    converged = False
    switch_count = 0
    for iteration in range(1, max_iter + 1):
        if log_mode:
            lv = logb - logsumexp(L + lu[:, None], axis=0)
            lu = loga - logsumexp(L + lv[None, :], axis=1)
        else:
            denom_v = K.T @ u
            if (denom_v <= 0).any() or not np.isfinite(denom_v).all():
                lu, lv = np.log(u), np.log(v)
                log_mode = True
                switch_count += 1
                continue
            v = b / denom_v
            denom_u = K @ v
            if (denom_u <= 0).any() or not np.isfinite(denom_u).all():
                lu, lv = np.log(u), np.log(v)
                log_mode = True
                switch_count += 1
                continue
            u = a / denom_u
            # The inverse rescaling leaves every P_ij exactly unchanged.
            if max(abs(np.log(u)).max(), abs(np.log(v)).max()) > 200:
                lu, lv = np.log(u), np.log(v)
                log_mode = True
                switch_count += 1
        if iteration % check_every == 0 or iteration == max_iter:
            if log_mode:
                r = np.exp(lu + logsumexp(L + lv[None, :], axis=1))
                c = np.exp(lv + logsumexp(L + lu[:, None], axis=0))
            else:
                r, c = u * (K @ v), v * (K.T @ u)
            max_abs = max(abs(r-a).max(), abs(c-b).max())
            max_rel = max((abs(r-a)/a).max(), (abs(c-b)/b).max())
            if max_abs <= abs_tol and max_rel <= rel_tol:
                converged = True
                break
    P = np.exp(L + lu[:, None] + lv[None, :]) if log_mode else u[:, None]*K*v[None, :]
    info = dict(iterations=iteration, converged=converged,
                stop_reason='marginal_tolerance' if converged else 'inner_iteration_cap',
                algorithm='log_domain' if log_mode else 'gauge_scaled_kernel',
                log_switches=switch_count, abs_tol=abs_tol, rel_tol=rel_tol,
                **plan_diagnostics(P, a, b))
    # Audit reconstructed P, not just the scaling-vector residual estimate.
    info['converged'] = bool(converged and info['finite'] and info['nonnegative']
        and max(info['row_max_abs'], info['col_max_abs']) <= abs_tol
        and max(info['row_max_relative'], info['col_max_relative']) <= rel_tol)
    if not info['converged'] and strict:
        raise TransportConvergenceError('Balanced Sinkhorn failed numerical acceptance', info)
    return P, info
