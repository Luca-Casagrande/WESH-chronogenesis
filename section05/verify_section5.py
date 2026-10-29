#!/usr/bin/env python3
"""Independent algebraic checks for WESH Section 5.

This verifies the declared finite realization and the explicitly added
product-law spatial sampling construction. It is not an empirical test,
not a derivation of the effective matching, and not a cosmological simulation.

Dependencies: numpy, sympy. Run: python section05/verify_section5.py
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import sympy as sp

q = sp.symbols('q', real=True)
nu, tau = sp.symbols('nu tau_T', positive=True)
gamma = sp.symbols('gamma_bi', nonnegative=True)
wl, wb = sp.symbols('w_loc w_bi', real=True)
scale = sp.symbols('scale', positive=True)
I = sp.eye(2)
Z = sp.diag(1, -1)
I4 = sp.eye(4)
Jx = (I4 + sp.kronecker_product(Z, I))/2
Jy = (I4 + sp.kronecker_product(I, Z))/2
L = Jx - Jy
jumps = [Jx, Jy, L]
coeffs = [nu, nu, gamma*q/2]
e = [I4[:, k] for k in range(4)]
vectors = [(e[0]+e[3])/sp.sqrt(2), (e[0]-e[3])/sp.sqrt(2),
           (e[1]+e[2])/sp.sqrt(2), (e[1]-e[2])/sp.sqrt(2)]
projectors = [v*v.T for v in vectors]
weights = sp.Matrix([(1+q)/4, (1+q)/4, (1-q)/4, (1-q)/4])
checks: list[str] = []

def equal(a, b, label: str) -> None:
    d = a-b
    if isinstance(d, sp.MatrixBase):
        ok = all(sp.simplify(z) == 0 for z in d)
    else:
        ok = sp.simplify(d) == 0
    if not ok:
        raise AssertionError(f'{label}: {sp.simplify(d)}')
    checks.append(label)

def dissipator(j: sp.Matrix, rho: sp.Matrix) -> sp.Matrix:
    return j*rho*j - (j*j*rho + rho*j*j)/2

rho = sum((weights[k]*projectors[k] for k in range(4)), sp.zeros(4))
equal(rho, (I4+q*sp.kronecker_product(Z,Z))/4, 'Bell barycenter')
equal(sum(weights), 1, 'Probability normalization')
rho_x = sp.Matrix(2, 2, lambda a,b: sum(rho[2*a+j,2*b+j] for j in range(2)))
rho_y = sp.Matrix(2, 2, lambda a,b: sum(rho[2*j+a,2*j+b] for j in range(2)))
equal(rho_x, I/2, 'Marginal x')
equal(rho_y, I/2, 'Marginal y')
# On q in [0,1], the difference from the product marginal has eigenvalues
# q/4, -q/4, -q/4, q/4, hence one half its trace norm is q/2.
equal(rho-sp.kronecker_product(rho_x,rho_y), q*sp.diag(1,-1,-1,1)/4,
      'Gate spectrum and product-faithfulness')
Q = sp.zeros(4)
state_rates: list[list[sp.Expr]] = []
for k,(v,pi) in enumerate(zip(vectors,projectors)):
    rates = []
    drift = sp.zeros(4)
    for a,(j,c) in enumerate(zip(jumps,coeffs)):
        mean = sp.trace(j*pi)
        delta = j-mean*I4
        variance = sp.simplify(sp.trace(pi*delta*delta))
        expected_var = sp.Rational(1,4) if a<2 else sp.Integer(k//2)
        equal(variance, expected_var, f'Variance state {k} channel {a}')
        equal(delta*delta*v, variance*v, f'No-event eigenidentity {k} {a}')
        drift += -c*(delta*delta*pi+pi*delta*delta)/2+c*variance*pi
        rates.append(c*variance)
        if variance != 0:
            update = delta*pi*delta/variance
            target = k^1
            equal(update, projectors[target], f'Centered update {k} {a}')
            equal(sp.trace(update*pi), 0, f'Orthogonality {k} {a}')
            Q[target,k] += c*variance
            Q[k,k] -= c*variance
    equal(drift, sp.zeros(4), f'Full no-event drift {k}')
    state_rates.append(rates)
equal(Q*weights, sp.zeros(4,1), 'Stationary full law')
equal(sp.ones(1,4)*Q, sp.zeros(1,4), 'Markov column sums')
for i in range(4):
    for j in range(i+1,4):
        equal(Q[i,j]*weights[j], Q[j,i]*weights[i], f'Detailed balance {i} {j}')
equal(sum((c*dissipator(j,rho) for c,j in zip(coeffs,jumps)),sp.zeros(4)),
      sp.zeros(4), 'Barycentric stationarity')
mean_bi = sum(weights[k]*state_rates[k][2] for k in range(4))
equal(mean_bi, gamma*q*(1-q)/4, 'Mean bilocal unit-lapse rate')
equal(tau*sum(weights[k]*sum(state_rates[k]) for k in range(4)),
      tau*(nu/2+gamma*q*(1-q)/4), 'Mean time-production')
response_weights = [wl,wl,wb]
proper_rates = [[sp.simplify(a*4/(tau*nu)) for a in row] for row in state_rates]
q2 = 0
joint_activity = 0
for k, row in enumerate(proper_rates):
    eps = k//2
    equal(row[0], 1/tau, f'Local proper rate x {k}')
    equal(row[1], 1/tau, f'Local proper rate y {k}')
    equal(row[2], 2*gamma*q*eps/(nu*tau), f'Bilocal proper rate {k}')
    beta = tau*sum(row)
    equal(beta, 2+2*gamma*q*eps/nu, f'Exact beta {k}')
    proportions = [a/sum(row) for a in row]
    W2 = sum(p*w*w for p,w in zip(proportions,response_weights))
    equal(beta*W2, 2*wl**2+2*gamma*q*wb**2*eps/nu, f'Joint beta-W2 {k}')
    q2 += weights[k]*sum(w*w*a for w,a in zip(response_weights,row))
    joint_activity += weights[k]*beta*W2
sigma2 = wl**2+gamma*wb**2*q*(1-q)/(2*nu)
equal(q2, 2*sigma2/tau, 'Weighted activity amplitude')
equal(joint_activity, 2*sigma2, 'Regional weighted-activity identity')
equal(sigma2.subs({nu:scale*nu,gamma:scale*gamma}),sigma2,
      'Common-rate cancellation at fixed law')
for endpoint in (0,1):
    equal(sigma2.subs(q,endpoint),wl**2,f'Local floor at q={endpoint}')
equal(sigma2.subs(q,sp.Rational(1,2)),wl**2+gamma*wb**2/(8*nu),
      'Maximum bilocal excess')
# The optional effective spatial realization is specified explicitly:
# disconnected copies, product initial law, independent drivers, identical q.
Q2 = sp.kronecker_product(Q,I4)+sp.kronecker_product(I4,Q)
equal(Q2*sp.kronecker_product(weights,weights),sp.zeros(16,1),
      'Stationarity of two-block product law')
p_odd = (1-q)/2
nb=5
Ef = sum(sp.binomial(nb,k)*p_odd**k*(1-p_odd)**(nb-k)*sp.Rational(k,nb)
         for k in range(nb+1))
Ef2 = sum(sp.binomial(nb,k)*p_odd**k*(1-p_odd)**(nb-k)*sp.Rational(k*k,nb*nb)
          for k in range(nb+1))
equal(Ef,p_odd,'Spatial fraction expectation (finite exact binomial check)')
equal(Ef2-Ef**2,(1-q*q)/(4*nb),'Spatial fraction variance')
equal(wl**2+gamma*q*wb**2*Ef/nu,sigma2,'Spatial average expectation')

# Independent numerical check of the finite constraint lift and its 2-state restriction.
rng=np.random.default_rng(20261014)
d=3
F=np.exp(2j*np.pi*np.outer(np.arange(d),np.arange(d))/d)/np.sqrt(d)
p=np.arange(-1,2,dtype=float)
T=F.conj().T@np.diag(p)@F
T2=T@T
Iq=np.eye(2)
C=np.kron(np.kron(np.diag(-p),Iq),np.eye(d))+np.kron(np.eye(d*2),np.diag(p))
V=np.zeros((d*2*d,d*2),complex)
for n in range(d):
    for a in range(2):
        V[(n*2+a)*d+n,n*2+a]=1
J=V@np.kron(T2,Iq)@V.conj().T
P=V@V.conj().T
Qcharge=np.kron(np.kron(np.eye(d),np.diag([-1.,1.])),np.eye(d))
U=F.conj().T@np.column_stack((np.array([1,0,1])/np.sqrt(2),np.array([0,1,0])))
encoded=V@np.kron(U,np.array([1,0]).reshape(2,1))
residuals={
    'constraint_isometry':float(np.linalg.norm(C@V)),
    'physical_support':float(np.linalg.norm(P@J@P-J)),
    'constraint_commutator':float(np.linalg.norm(C@J-J@C)),
    'neutral_charge_commutator':float(np.linalg.norm(Qcharge@J-J@Qcharge)),
    'invariant_two_state_subspace':float(np.linalg.norm(J@encoded-encoded@np.diag([1.,0.]))),
    'two_state_J':float(np.linalg.norm(encoded.conj().T@J@encoded-np.diag([1.,0.]))),
}
for name,res in residuals.items():
    if res>1e-12:
        raise AssertionError(f'{name}: {res}')
# Independent direct 4x4 state calculations for random parameter choices.
Jn=[np.array(j,dtype=complex) for j in jumps]
pin=[np.array(pi,dtype=complex) for pi in projectors]
max_stationarity=0.
max_amplitude=0.
for _ in range(200):
    qq=rng.uniform(0,1); nn=rng.uniform(.1,5); gg=rng.uniform(0,7)
    aa,bb=rng.normal(size=2); tt=rng.uniform(.1,4)
    probs=np.array([(1+qq)/4]*2+[(1-qq)/4]*2)
    rn=sum(pr*pi for pr,pi in zip(probs,pin))
    cc=[nn,nn,gg*qq/2]
    stat=sum(c*(j@rn@j-(j@j@rn+rn@j@j)/2) for c,j in zip(cc,Jn))
    max_stationarity=max(max_stationarity,float(np.linalg.norm(stat)))
    actual=0.
    for pr,pi in zip(probs,pin):
        for j,c,w in zip(Jn,cc,[aa,aa,bb]):
            mean=np.trace(pi@j).real
            var=(np.trace(pi@j@j).real-mean*mean)
            actual += pr*w*w*c*var*4/(tt*nn)
    pred=aa*aa+gg/nn*bb*bb*qq*(1-qq)/2
    max_amplitude=max(max_amplitude,abs(actual*tt/2-pred))
    assert aa*aa-1e-12 <= pred <= aa*aa+gg/nn*bb*bb/8+1e-12
assert max_stationarity<1e-12 and max_amplitude<1e-12
report={
 'all_checks_passed':True,
 'exact_symbolic_checks':len(checks),
 'checks':checks,
 'finite_lift_residuals':residuals,
 'randomized_direct_checks':200,
 'maximum_stationarity_residual':max_stationarity,
 'maximum_amplitude_residual':max_amplitude,
 'scope':'Algebraic verification of the stated finite model. Product-law spatial sampling is an explicitly declared sufficient effective assignment. No empirical or cosmological data are used.'
}
output=Path(__file__).with_name('VERIFICHE_ALGEBRICHE.json')
output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))
