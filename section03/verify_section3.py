#!/usr/bin/env python3
"""Finite-dimensional regression checks for WESH Section 3.

These checks audit algebraic identities and limiting conventions. They are
neither proofs in place of the manuscript's arguments nor empirical tests
or simulations of a complete cosmology. No network access is required.
Dependencies: numpy, scipy, sympy. Run: python verify_section3.py
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm
import sympy as sp

TOL = 2e-10
RNG = np.random.default_rng(31415926)
RES: dict[str, float | int | bool | list[float] | str] = {}


def adj(a: np.ndarray) -> np.ndarray:
    return a.conj().T


def herm(a: np.ndarray) -> np.ndarray:
    return (a + adj(a)) / 2


def tnorm(a: np.ndarray) -> float:
    return float(np.linalg.svd(a, compute_uv=False).sum())


def density(d: int) -> np.ndarray:
    x = RNG.normal(size=(d, d)) + 1j * RNG.normal(size=(d, d))
    r = x @ adj(x)
    return r / np.trace(r)


def pure(d: int) -> tuple[np.ndarray, np.ndarray]:
    p = RNG.normal(size=d) + 1j * RNG.normal(size=d)
    p = p / np.linalg.norm(p)
    return p, np.outer(p, p.conj())


def reduce(a: np.ndarray, dims: Sequence[int], keep: Sequence[int]) -> np.ndarray:
    """Partial trace, retaining factors in their original order."""
    if len(set(keep)) != len(keep) or any(k < 0 or k >= len(dims) for k in keep):
        raise ValueError('invalid subsystem indices')
    shape = list(dims)
    if a.shape != (int(np.prod(shape)), int(np.prod(shape))):
        raise ValueError('matrix and factor dimensions disagree')
    t = a.reshape(shape + shape)
    for k in sorted(set(range(len(shape))) - set(keep), reverse=True):
        t = np.trace(t, axis1=k, axis2=k + len(shape))
        shape.pop(k)
    d = int(np.prod(shape))
    return t.reshape(d, d)


def gate(r: np.ndarray, dims: Sequence[int], x: int, y: int) -> float:
    x, y = sorted((x, y))
    pair = reduce(r, dims, (x, y))
    return 0.5 * tnorm(pair - np.kron(reduce(r, dims, (x,)), reduce(r, dims, (y,))))


def dissip(j: np.ndarray, r: np.ndarray) -> np.ndarray:
    j2 = j @ j
    return j @ r @ j - (j2 @ r + r @ j2) / 2


def gen(r: np.ndarray, js: Sequence[np.ndarray], cs: Sequence[float], h: np.ndarray) -> np.ndarray:
    out = -1j * (h @ r - r @ h)  # hbar=1 only in these dimensionless checks
    for j, c in zip(js, cs):
        out += c * dissip(j, r)
    return out


def event_data(pi: np.ndarray, js: Sequence[np.ndarray], cs: Sequence[float], h: np.ndarray):
    ident = np.eye(pi.shape[0], dtype=complex)
    drift = -1j * (h @ pi - pi @ h)
    vector_coeff = -1j * h.copy()
    updates, rates, variances = [], [], []
    for j, c in zip(js, cs):
        mean = float(np.trace(pi @ j).real)
        d = j - mean * ident
        v = max(0.0, float(np.trace(pi @ d @ d).real))
        u = d @ pi @ d / v if v > 1e-14 else pi.copy()
        drift += c * (-(d @ d @ pi + pi @ d @ d) / 2 + v * pi)
        vector_coeff -= c * (d @ d - v * ident) / 2
        updates.append(u)
        rates.append(c * v)
        variances.append(v)
    return drift, vector_coeff, updates, np.array(rates), np.array(variances)


def check(name: str, value: float, tol: float = TOL) -> None:
    RES[name] = float(value)
    if not np.isfinite(value) or value > tol:
        raise AssertionError(f'{name}: {value} > {tol}')


I2 = np.eye(2, dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
J = (I2 + Z) / 2
j0, j1 = np.kron(J, I2), np.kron(I2, J)
JPAIR = [j0, j1, j0 - j1]
NU, GAMMA = 1.3, 0.7

# 1. Trace-norm gate: product zero, bounds, Lipschitz, unitary invariance,
#    and contraction under a product channel.
max_lip = max_unitary = max_contract = max_product = 0.0
for _ in range(40):
    r, e = density(4), density(4)
    cr, ce = gate(r, (2, 2), 0, 1), gate(e, (2, 2), 0, 1)
    max_lip = max(max_lip, abs(cr - ce) - 1.5 * tnorm(r - e))
    ra, rb = density(2), density(2)
    max_product = max(max_product, gate(np.kron(ra, rb), (2, 2), 0, 1))
    u1, _ = np.linalg.qr(RNG.normal(size=(2, 2)) + 1j * RNG.normal(size=(2, 2)))
    u2, _ = np.linalg.qr(RNG.normal(size=(2, 2)) + 1j * RNG.normal(size=(2, 2)))
    u = np.kron(u1, u2)
    max_unitary = max(max_unitary, abs(gate(u @ r @ adj(u), (2, 2), 0, 1) - cr))
    p = 0.41
    ks = [np.diag([1, np.sqrt(1-p)]), np.array([[0, np.sqrt(p)], [0, 0]])]
    out = sum(np.kron(k, I2) @ r @ adj(np.kron(k, I2)) for k in ks)
    max_contract = max(max_contract, gate(out, (2, 2), 0, 1) - cr)
    assert -TOL <= cr <= 1 + TOL
check('gate_lipschitz_violation', max_lip)
check('gate_product_zero', max_product)
check('gate_local_unitary_residual', max_unitary)
check('gate_local_channel_contraction_violation', max_contract)

# 2. Centered update, norm preservation, and cancellation of conditional
#    drift with mean jumps. No ensemble simulation is used.
mean_res = norm_res = orth_res = scale_res = 0.0
for d in (2, 3, 5):
    for _ in range(20):
        psi, pi = pure(d)
        js = [herm(RNG.normal(size=(d,d))+1j*RNG.normal(size=(d,d))) for _ in range(3)]
        h = herm(RNG.normal(size=(d,d))+1j*RNG.normal(size=(d,d)))
        cs = [0.2, 1.1, 0.0]
        b, bv, us, rates, vs = event_data(pi, js, cs, h)
        calculated = b + sum(l * (u - pi) for l, u in zip(rates, us))
        mean_res = max(mean_res, tnorm(calculated - gen(pi, js, cs, h)))
        norm_res = max(norm_res, abs(float(2*np.vdot(psi, bv @ psi).real)))
        for v, u in zip(vs, us):
            if v > 1e-14:
                orth_res = max(orth_res, abs(np.trace(pi @ u)), abs(np.trace(u)-1))
        scaled_js = [-2.3*j + 0.71*np.eye(d) for j in js]
        bs, _, uss, ls, _ = event_data(pi, scaled_js, [c/2.3**2 for c in cs], h)
        scale_res = max(scale_res, tnorm(bs-b), np.max(abs(ls-rates)),
                        max(tnorm(u-v) for u,v in zip(us,uss)))
check('conditional_drift_identity', mean_res)
check('no_event_norm_derivative', norm_res)
check('orthogonality_and_normalization', float(orth_res))
check('identity_shift_and_rate_rescaling', scale_res)

# 3. Full covariance, including sign reversal of the edge representative
#    when the two endpoints are interchanged.
swap = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], complex)
h = 0.3*(np.kron(X,I2)+np.kron(I2,X)) + 0.1*np.kron(Z,Z)
rho = density(4)
_, pi = pure(4)
cs = [NU, NU, GAMMA*gate(rho,(2,2),0,1)]
b, _, us, ls, _ = event_data(pi, JPAIR, cs, h)
rhog, pig = swap@rho@swap, swap@pi@swap
csg = [NU,NU,GAMMA*gate(rhog,(2,2),0,1)]
bg, _, ugs, lgs, _ = event_data(pig,JPAIR,csg,h)
check('covariant_no_event_drift', tnorm(bg-swap@b@swap))
check('covariant_channel_intensities', float(np.max(abs(lgs-ls[[1,0,2]]))))
check('covariant_post_event_projectors', max(tnorm(ugs[a]-swap@us[q]@swap)
                                           for a,q in enumerate((1,0,2))))

# 4. Stationarity on the Bell-sign law and the distinction between quantum
#    and ensemble variance.
bell = [np.array(v, complex)/np.sqrt(2) for v in ((1,0,0,1),(1,0,0,-1),(0,1,1,0),(0,1,-1,0))]
pis = [np.outer(v,v.conj()) for v in bell]
q = 0.37
weights = np.array([(1+q)/4]*2+[(1-q)/4]*2)
rho = sum(w*p for w,p in zip(weights,pis))
cs = [NU, NU, GAMMA*q/2]
rate_matrix = np.zeros((4,4))
activity = []
for k, pi in enumerate(pis):
    b, _, us, ls, vs = event_data(pi,JPAIR,cs,np.zeros((4,4),complex))
    check(f'bell_no_event_drift_{k}', tnorm(b))
    for l, u in zip(ls,us):
        if l > 1e-14:
            target = int(np.argmin([tnorm(u-p) for p in pis]))
            assert tnorm(u-pis[target]) < TOL
            rate_matrix[target,k] += l
            rate_matrix[k,k] -= l
    activity.append(float(ls.sum()))
check('bell_law_stationarity', float(np.linalg.norm(rate_matrix@weights)))
check('bell_barycenter_stationarity', tnorm(gen(rho,JPAIR,cs,np.zeros((4,4),complex))))
check('bell_expected_total_activity', abs(float(weights@np.array(activity))-(NU/2+GAMMA*q*(1-q)/4)))
RES['stationary_path_intensities_even_odd'] = [activity[0],activity[2]]
RES['stationary_expected_intensity'] = float(weights@np.array(activity))
for a,j in enumerate(JPAIR):
    jbar = float(np.trace(rho@j).real)
    vlaw = 0.; vcond = 0.
    for w,pi in zip(weights,pis):
        jm=float(np.trace(pi@j).real)
        vcond += w*float((np.trace(pi@j@j)-jm*jm).real)
        vlaw += w*(jm-jbar)**2
    check(f'law_quantum_variance_split_{a}', abs(float(np.trace(rho@j@j).real)-jbar*jbar-vlaw-vcond))
# A second law with the same barycenter consists of computational eigenstates.
for k in range(4):
    v=np.eye(4,dtype=complex)[:,k]; pi=np.outer(v,v.conj())
    b,_,_,ls,_=event_data(pi,JPAIR,cs,np.zeros((4,4),complex))
    check(f'sharp_null_activity_{k}', float(ls.sum()))
    check(f'sharp_null_drift_{k}', tnorm(b))
check('same_barycenter_computational_law', tnorm(rho-np.diag(np.diag(rho))))

# 5. Spectator extension and frozen CPTP propagation.
joint = density(8)
sys = reduce(joint,(4,2),(0,))
h = np.kron(X,I2)+0.2*np.kron(I2,X)
cs = [NU,NU,GAMMA*gate(sys,(2,2),0,1)]
jjs=[np.kron(j,I2) for j in JPAIR]
djoint=gen(joint,jjs,cs,np.kron(h,I2))
check('ancilla_system_marginal',tnorm(reduce(djoint,(4,2),(0,))-gen(sys,JPAIR,cs,h)))
check('ancilla_spectator_marginal',tnorm(reduce(djoint,(4,2),(1,))))
d=4; ident=np.eye(d)
mat=-1j*(np.kron(ident,h)-np.kron(h.T,ident))
for j,c in zip(JPAIR,cs):
    jj=j@j
    mat += c*(np.kron(j.conj(),j)-0.5*np.kron(ident,jj)-0.5*np.kron(jj.T,ident))
out=(expm(0.41*mat)@sys.reshape(-1,order='F')).reshape((d,d),order='F')
check('frozen_trace',abs(np.trace(out)-1))
check('frozen_hermiticity',tnorm(out-adj(out)))
check('frozen_positivity_violation',max(0.,-float(np.linalg.eigvalsh(herm(out)).min())))

# 6. Reduced instantaneous dynamics on a three-vertex path, before/after a
#    nonselective operation outside A^+ (A={0}, A^+={0,1}, B={2}).
loc=[]
for k in range(3):
    ops=[I2,I2,I2]; ops[k]=J
    loc.append(np.kron(np.kron(ops[0],ops[1]),ops[2]))
js=loc+[loc[0]-loc[1],loc[1]-loc[2]]
h=0.17*(np.kron(np.kron(X,X),I2)+np.kron(np.kron(I2,X),X))
rho=density(8)
ks=[np.diag([1,np.sqrt(0.6)]),np.array([[0,np.sqrt(0.4)],[0,0]])]
rhop=sum(np.kron(np.eye(4),k)@rho@adj(np.kron(np.eye(4),k)) for k in ks)
def chain_coeff(r):
    return [NU]*3+[GAMMA*0.8*gate(r,(2,2,2),0,1),GAMMA*1.2*gate(r,(2,2,2),1,2)]
a=reduce(gen(rho,js,chain_coeff(rho),h),(2,2,2),(0,))
b=reduce(gen(rhop,js,chain_coeff(rhop),h),(2,2,2),(0,))
check('graph_local_instantaneous_influence',tnorm(a-b))

# 7. Reparametrization: integrated intensity and its transformed image.
S=1.7
integrand=lambda s:(1+s)*(1+0.2*np.sin(s))
first=quad(integrand,0,S,epsabs=1e-12)[0]
# s'=exp(s)-1, ds/ds'=1/(1+s')
second=quad(lambda z:integrand(np.log1p(z))/(1+z),0,np.expm1(S),epsabs=1e-12)[0]
check('compensator_reparametrization',abs(first-second))

# 8. Exact scalar ingredient of the martingale strong-law proof.
a = sp.symbols('a', nonnegative=True)
expr = sp.integrate(1/(1+a)**2,(a,0,a))
assert sp.simplify(expr-(1-1/(1+a))) == 0
RES['martingale_weight_integral_symbolic'] = True

# 9. Centered jumps on the finite constraint sector from Section 2.
d=3; qdim=2
p=np.array([-1.,0.,1.]); F=np.exp(2j*np.pi*np.outer(np.arange(d),np.arange(d))/d)/np.sqrt(d)
t=adj(F)@np.diag(p)@F; t2=herm(t@t)
C=np.kron(np.kron(np.diag(-p),np.eye(qdim)),np.eye(d))+np.kron(np.eye(d*qdim),np.diag(p))
V=np.zeros((d*qdim*d,d*qdim),complex)
for n in range(d):
    for k in range(qdim): V[(n*qdim+k)*d+n,n*qdim+k]=1
P=V@adj(V)
Jlift=V@np.kron(t2,np.eye(qdim))@adj(V)
charge=np.kron(np.kron(np.eye(d),Z),np.eye(d))
psi,pi0=pure(d*qdim); pi=V@pi0@adj(V)
b,_,us,ls,_=event_data(pi,[Jlift],[1.2],np.zeros_like(C))
check('physical_centered_drift_support',tnorm((np.eye(len(C))-P)@b))
check('physical_event_support',tnorm((np.eye(len(C))-P)@us[0]))
check('noether_mean_drift',abs(np.trace(charge@gen(pi,[Jlift],[1.2],np.zeros_like(C)))))

RES['status']='PASS'
RES['seed']=31415926
RES['scope']='Algebraic and numerical regression checks; no empirical or full-cosmological simulation.'
outpath=Path(__file__).with_name('section3_checks.json')
outpath.write_text(json.dumps(RES,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print(json.dumps(RES,indent=2,ensure_ascii=False))
