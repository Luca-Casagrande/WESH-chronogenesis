#!/usr/bin/env python3
"""Deterministic algebraic checks for WESH Section 2.

Checks the equations written in the section. It is not a simulation of
chronogenesis or a test of continuum Wheeler--DeWitt dynamics. The spectral
construction is restated from the project's physical-realization report;
this script is an independent implementation of those formulas.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

TOL = 3e-11
RNG = np.random.default_rng(20261013)
RESULTS: dict[str, float | int | str] = {}


def adj(a: np.ndarray) -> np.ndarray:
    return a.conj().T


def comm(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a @ b - b @ a


def diss(j: np.ndarray, r: np.ndarray) -> np.ndarray:
    return j @ r @ j - (j @ j @ r + r @ j @ j) / 2


def herm(d: int) -> np.ndarray:
    a = RNG.normal(size=(d, d)) + 1j * RNG.normal(size=(d, d))
    return (a + adj(a)) / 2


def check(name: str, residual: np.ndarray | float | complex) -> None:
    v = float(np.linalg.norm(residual))
    RESULTS[name] = v
    if not np.isfinite(v) or v > TOL:
        raise AssertionError(f'{name}: residual {v:.5g} exceeds {TOL}')


def sqrt_psd(a: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eigh((a + adj(a)) / 2)
    if float(vals.min()) < -TOL:
        raise ValueError('Nonpositive Kraus complement')
    return (vecs * np.sqrt(np.maximum(vals, 0))) @ adj(vecs)


def noether_checks() -> None:
    for i in range(5):
        d = 5
        h, q, j1, j2 = (herm(d) for _ in range(4))
        rates = (0.7, 1.3)
        ladj = 1j * comm(h, q) + sum(c * diss(j, q) for c,j in zip(rates,(j1,j2)))
        expected = -0.5 * sum(c * np.linalg.norm(comm(j,q))**2 for c,j in zip(rates,(j1,j2)))
        check(f'noether_trace_identity_{i}', np.trace(q @ ladj)-expected)
        check(f'double_commutator_{i}', diss(j1,q)+comm(j1,comm(j1,q))/2)
        b = herm(d)
        check(f'trace_adjoint_sign_{i}', np.trace(q@(-1j*comm(h,b)))-np.trace((1j*comm(h,q))@b))
        check(f'real_identity_shift_{i}', diss(j1+1.37*np.eye(d),b)-diss(j1,b))


def spectral_checks(d: int) -> None:
    qdim = 2
    p = np.arange(-(d//2), d//2+1,dtype=float)
    nodes = p/(d//2)
    f = np.exp(2j*np.pi*np.outer(np.arange(d),np.arange(d))/d)/np.sqrt(d)
    t = adj(f) @ np.diag(nodes) @ f
    j = t @ t
    iq = np.eye(qdim)
    it = np.eye(d)
    qvals = np.diag([-1.,1.])
    c = np.kron(np.kron(np.diag(-p),iq),it)+np.kron(np.eye(d*qdim),np.diag(p))
    v = np.zeros((d*qdim*d,d*qdim),complex)
    for n in range(d):
        for q in range(qdim):
            v[(n*qdim+q)*d+n,n*qdim+q]=1.
    proj = v@adj(v)
    jr = np.kron(j,iq)
    jl = v@jr@adj(v)
    qkin = np.kron(np.kron(it,qvals),it)
    bare = np.kron(np.eye(d*qdim),j)
    mu2 = np.mean(nodes**2)
    prefix = f'd{d}_'
    check(prefix+'isometry',adj(v)@v-np.eye(d*qdim))
    check(prefix+'constraint_embedding',c@v)
    check(prefix+'constraint_commutator',comm(c,jl))
    check(prefix+'noether_commutator',comm(qkin,jl))
    check(prefix+'intertwining',jl@v-v@jr)
    check(prefix+'bare_compression_central',adj(v)@bare@v-mu2*np.eye(d*qdim))
    a,b = herm(d),herm(d)
    lift=lambda z: v@np.kron(z,iq)@adj(v)
    check(prefix+'homomorphism',lift(a@b)-lift(a)@lift(b))
    check(prefix+'star_homomorphism',lift(adj(a))-adj(lift(a)))
    check(prefix+'physical_unit',lift(it)-proj)
    vals,vecs=np.linalg.eigh(jr)
    od=np.outer(vecs[:,0],vecs[:,-1].conj())
    check(prefix+'nontrivial_dissipator_identity',diss(jr,od)+0.5*(vals[0]-vals[-1])**2*od)
    RESULTS[prefix+'spectral_gap'] = float(vals[-1]-vals[0])
    if vals[-1]-vals[0] < 0.1: raise AssertionError('Trivial jump')
    psi=RNG.normal(size=d*qdim)+1j*RNG.normal(size=d*qdim)
    psi/=np.linalg.norm(psi)
    rho=np.outer(psi,psi.conj())
    delta=jr-np.trace(jr@rho).real*np.eye(d*qdim)
    bdelta=0.8*(delta@delta)
    eps=0.2/max(np.linalg.norm(bdelta,2),1.)
    k1=np.sqrt(eps*0.8)*delta
    k0=sqrt_psd(np.eye(d*qdim)-eps*bdelta)
    check(prefix+'kraus_completeness',adj(k0)@k0+adj(k1)@k1-np.eye(d*qdim))
    qr=np.kron(it,qvals)
    check(prefix+'kraus_noether',adj(k0)@qr@k0+adj(k1)@qr@k1-qr)
    for n,k in enumerate((k0,k1)):
        kl=v@k@adj(v)
        check(prefix+f'kraus_branch_closure_{n}',(np.eye(len(c))-proj)@kl@proj)
    check(prefix+'centered_generator_identity',diss(delta,rho)-diss(jr,rho))
    if d==3:
        even=np.array([1.,0.,1.])/np.sqrt(2)
        zero=np.array([0.,1.,0.])
        w=adj(f)@np.column_stack([even,zero])
        reversal=adj(f)@np.fliplr(np.eye(d))@f
        jqubit=np.diag([1.,0.])
        check('qubit_restriction_isometry',adj(w)@w-np.eye(2))
        check('qubit_restriction_invariance',j@w-w@jqubit)
        check('qubit_orientation_neutrality',reversal@w-w)
        check('qubit_projector_form',adj(w)@j@w-jqubit)
        jx=np.kron(jqubit,np.eye(2)); jy=np.kron(np.eye(2),jqubit)
        check('two_factor_local_commutation',comm(jx,jy))
        swap=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1.]])
        check('bilocal_endpoint_sign',swap@(jx-jy)@swap+(jx-jy))


def main() -> None:
    noether_checks()
    for d in (3,5,7):
        spectral_checks(d)
    residuals=[v for k,v in RESULTS.items() if isinstance(v,float) and not k.endswith('spectral_gap')]
    RESULTS['maximum_checked_residual']=max(residuals)
    RESULTS['status']='PASS_ALGEBRAIC_CHECKS'
    RESULTS['scope']='Finite algebra only; no continuum or empirical validation.'
    dest=Path(__file__).with_name('section2_checks.json')
    dest.write_text(json.dumps(RESULTS,indent=2)+'\n')
    print(json.dumps({k:RESULTS[k] for k in ('status','maximum_checked_residual','scope')},indent=2))

if __name__=='__main__':
    main()
