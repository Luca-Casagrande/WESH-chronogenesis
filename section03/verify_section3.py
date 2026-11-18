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


def negativity(pair: np.ndarray, dx: int, dy: int) -> float:
    """Negativity of a pair state, capped at one: min(1, ||pair^{T_x}||_1 - 1)."""
    pt = pair.reshape(dx, dy, dx, dy).transpose(2, 1, 0, 3).reshape(dx*dy, dx*dy)
    return min(1.0, tnorm(pt) - 1.0)


def upsilon(pi: np.ndarray, dims: Sequence[int], x: int, y: int) -> np.ndarray:
    """Negativity-weighted pair correlation of one conditional state."""
    x, y = sorted((x, y))
    pair = reduce(pi, dims, (x, y))
    product = np.kron(reduce(pi, dims, (x,)), reduce(pi, dims, (y,)))
    return negativity(pair, dims[x], dims[y]) * (pair - product)


def gate(law, dims: Sequence[int], x: int, y: int) -> float:
    """Entanglement gate of a law given as (weight, pure state) pairs."""
    return 0.5 * tnorm(sum(w * upsilon(pi, dims, x, y) for w, pi in law))


def product_form_gate(law, dims: Sequence[int], x: int, y: int) -> float:
    """Mean negativity times barycentric correlation; comparison only."""
    x, y = sorted((x, y))
    bary = sum(w * pi for w, pi in law)
    pair = reduce(bary, dims, (x, y))
    corr = 0.5 * tnorm(pair - np.kron(reduce(bary, dims, (x,)), reduce(bary, dims, (y,))))
    mean = sum(w * negativity(reduce(pi, dims, (x, y)), dims[x], dims[y]) for w, pi in law)
    return mean * corr


def random_law(d: int, count: int):
    weights = RNG.dirichlet(np.ones(count))
    return [(float(w), pure(d)[1]) for w in weights]


def ket(*amplitudes) -> np.ndarray:
    v = np.array(amplitudes, dtype=complex)
    v = v / np.linalg.norm(v)
    return np.outer(v, v.conj())


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

# 1. Entanglement gate of a law: range, zero set, cancellation, Lipschitz
#    bound in the coupled (Wasserstein) sense, symmetry and local unitaries.
SWAP2 = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], complex)
BELLS = [ket(1,0,0,1), ket(1,0,0,-1), ket(0,1,1,0), ket(0,1,-1,0)]
max_range = max_point = max_unitary = max_swap = 0.0
max_state_lip = max_law_lip = 0.0
for _ in range(40):
    law = random_law(4, 3)
    c = gate(law, (2, 2), 0, 1)
    max_range = max(max_range, -c, c - 1)
    _, pi = pure(4)
    eta = negativity(pi, 2, 2)
    # On a point law the gate is positive exactly when the negativity is.
    max_point = max(max_point, abs(gate([(1.0, pi)], (2, 2), 0, 1)
                                   - eta * 0.5 * tnorm(pi - np.kron(reduce(pi,(2,2),(0,)), reduce(pi,(2,2),(1,))))))
    u1, _ = np.linalg.qr(RNG.normal(size=(2, 2)) + 1j * RNG.normal(size=(2, 2)))
    u2, _ = np.linalg.qr(RNG.normal(size=(2, 2)) + 1j * RNG.normal(size=(2, 2)))
    u = np.kron(u1, u2)
    max_unitary = max(max_unitary, abs(gate([(w, u @ p @ adj(u)) for w, p in law], (2, 2), 0, 1) - c))
    max_swap = max(max_swap, abs(gate([(w, SWAP2 @ p @ SWAP2) for w, p in law], (2, 2), 0, 1) - c))
    # Per-state Lipschitz constant d_xy + 3/2, also inside a larger system.
    for dims, pair in (((2, 2), (0, 1)), ((3, 2), (0, 1)), ((2, 2, 2), (0, 2))):
        d = int(np.prod(dims)); dxy = min(dims[pair[0]], dims[pair[1]])
        _, p1 = pure(d)
        v = RNG.normal(size=d) + 1j * RNG.normal(size=d)
        w0 = np.linalg.eigh(p1)[1][:, -1] + RNG.choice([1e-2, 1e-1, 1.0]) * v
        p2 = np.outer(w0, w0.conj()) / np.vdot(w0, w0).real
        lhs = 0.5 * tnorm(upsilon(p1, dims, *pair) - upsilon(p2, dims, *pair))
        max_state_lip = max(max_state_lip, lhs - (dxy + 1.5) * tnorm(p1 - p2))
    # Two laws coupled state by state: W_1 is at most the mean pair distance.
    other = []
    for w, p1 in law:
        v = RNG.normal(size=4) + 1j * RNG.normal(size=4)
        w0 = np.linalg.eigh(p1)[1][:, -1] + 0.3 * v
        other.append((w, np.outer(w0, w0.conj()) / np.vdot(w0, w0).real))
    coupling_cost = sum(w * tnorm(p1 - p2) for (w, p1), (_, p2) in zip(law, other))
    max_law_lip = max(max_law_lip, abs(c - gate(other, (2, 2), 0, 1)) - 3.5 * coupling_cost)
check('gate_range_violation', max(0.0, max_range))
check('gate_point_law_identity', max_point)
check('gate_local_unitary_residual', max_unitary)
check('gate_exchange_residual', max_swap)
check('gate_state_lipschitz_violation', max(0.0, max_state_lip))
check('gate_law_lipschitz_violation', max(0.0, max_law_lip))
# Unentangled realizations drop out before the average is taken.
plus, minus = np.array([1, 1]) / np.sqrt(2), np.array([1, -1]) / np.sqrt(2)
classical = [(0.5, ket(*np.kron(plus, plus))), (0.5, ket(*np.kron(minus, minus)))]
check('gate_classical_mixture_zero', gate(classical, (2, 2), 0, 1))
RES['classical_mixture_barycentric_correlation'] = float(
    0.5 * tnorm(sum(w * p for w, p in classical) - np.eye(4) / 4))
# Unpolarized Bell realizations plus classically correlated product ones.
mixed = [(0.125, b) for b in BELLS] + [(0.25, ket(1,0,0,0)), (0.25, ket(0,0,0,1))]
check('gate_mixed_branches_zero', gate(mixed, (2, 2), 0, 1))
check('product_form_comparison_value', abs(product_form_gate(mixed, (2, 2), 0, 1) - 0.125))
# A pair whose reduction is separable because each site is entangled elsewhere.
ghz = ket(1,0,0,0,0,0,0,1)
check('gate_third_party_separable_zero', gate([(1.0, ghz)], (2, 2, 2), 0, 1))
# Bell law: unit negativity, gate q/2.
for qq in (0.0, 0.37, 1.0):
    bell_law = list(zip([(1+qq)/4]*2 + [(1-qq)/4]*2, BELLS))
    check(f'gate_bell_law_q={qq}', abs(gate(bell_law, (2, 2), 0, 1) - qq / 2))
check('bell_unit_negativity', max(abs(negativity(b, 2, 2) - 1) for b in BELLS))

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
law = random_law(4, 3)
_, pi = pure(4)
cs = [NU, NU, GAMMA*gate(law,(2,2),0,1)]
b, _, us, ls, _ = event_data(pi, JPAIR, cs, h)
lawg, pig = [(w, swap@p@swap) for w, p in law], swap@pi@swap
csg = [NU,NU,GAMMA*gate(lawg,(2,2),0,1)]
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
cs = [NU, NU, GAMMA*gate(list(zip(weights,pis)),(2,2),0,1)]
check('bell_gate_value', abs(cs[2]-GAMMA*q/2))
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

# 5. Spectator extension and frozen CPTP propagation. The conditional
#    states are pure states of system (two qubits) and spectator (one qubit);
#    the gate reads only their pair reductions on the system.
joint_law = random_law(8, 3)
joint = sum(w*p for w, p in joint_law)
sys = reduce(joint,(4,2),(0,))
h = np.kron(X,I2)+0.2*np.kron(I2,X)
cgate = gate(joint_law,(2,2,2),0,1)
cs = [NU,NU,GAMMA*cgate]
jjs=[np.kron(j,I2) for j in JPAIR]
djoint=gen(joint,jjs,cs,np.kron(h,I2))
check('ancilla_system_marginal',tnorm(reduce(djoint,(4,2),(0,))-gen(sys,JPAIR,cs,h)))
check('ancilla_spectator_marginal',tnorm(reduce(djoint,(4,2),(1,))))
# A unitary on the spectator leaves every conditional pair reduction unchanged.
ur, _ = np.linalg.qr(RNG.normal(size=(2, 2)) + 1j * RNG.normal(size=(2, 2)))
ur = np.kron(np.eye(4), ur)
check('ancilla_isometry_gate_invariance',
      abs(gate([(w, ur@p@adj(ur)) for w, p in joint_law],(2,2,2),0,1)-cgate))
# Resolving the spectator into recorded outcomes is a different law.
purified = np.zeros(8, complex)
purified[[0, 6]] = 1 / 2            # |00>|0>, |11>|0>
purified[[3, 5]] = 1 / 2            # |01>|1>, |10>|1>
purified = np.outer(purified, purified.conj())
check('unresolved_spectator_gate', gate([(1.0, purified)],(2,2,2),0,1))
resolved = [(0.5, BELLS[0]), (0.5, BELLS[2])]
check('resolved_spectator_gate', abs(gate(resolved,(2,2),0,1)-0.5))
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
#    unitary operation outside A^+ (A={0}, A^+={0,1}, B={2}), applied to
#    every conditional state of the law.
loc=[]
for k in range(3):
    ops=[I2,I2,I2]; ops[k]=J
    loc.append(np.kron(np.kron(ops[0],ops[1]),ops[2]))
js=loc+[loc[0]-loc[1],loc[1]-loc[2]]
h=0.17*(np.kron(np.kron(X,X),I2)+np.kron(np.kron(I2,X),X))
law3=random_law(8,3)
ub,_=np.linalg.qr(RNG.normal(size=(2,2))+1j*RNG.normal(size=(2,2)))
ub=np.kron(np.eye(4),ub)
law3p=[(w,ub@p@adj(ub)) for w,p in law3]
def chain_coeff(law):
    return [NU]*3+[GAMMA*0.8*gate(law,(2,2,2),0,1),GAMMA*1.2*gate(law,(2,2,2),1,2)]
def bary(law):
    return sum(w*p for w,p in law)
a=reduce(gen(bary(law3),js,chain_coeff(law3),h),(2,2,2),(0,))
b=reduce(gen(bary(law3p),js,chain_coeff(law3p),h),(2,2,2),(0,))
check('graph_local_instantaneous_influence',tnorm(a-b))
check('graph_local_edge_gate_unchanged',
      abs(gate(law3,(2,2,2),0,1)-gate(law3p,(2,2,2),0,1)))

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
