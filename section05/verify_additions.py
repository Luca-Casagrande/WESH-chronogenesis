#!/usr/bin/env python3
"""Exact checks for the Section 5 regional identity and q=1 endpoint.

No observations are used. The volume assignment, full collection and the
product initial law are hypotheses of the effective realization, not outputs
of these algebraic checks.
"""
from pathlib import Path
import importlib.util
import json
import sympy as sp

here=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('base_checks',here/'verify_section5.py')
base=importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
passed=[]
def check(a,b,name):
    diff=a-b
    terms=list(diff) if isinstance(diff,sp.MatrixBase) else [diff]
    if not all(sp.simplify(x)==0 for x in terms):
        raise AssertionError(name)
    passed.append(name)

q,r=sp.symbols('q r_bi', nonnegative=True)
n=sp.symbols('n_b',integer=True,positive=True)
LP,dt,nu=sp.symbols('L_P Delta_tau nu',positive=True)
wl,wb=sp.symbols('w_loc w_bi',real=True)
K=sp.symbols('K',integer=True,nonnegative=True)
sig=wl**2+r*wb**2*q*(1-q)/2
V=2*n*LP**3*dt
rate=2*(n*wl**2+r*q*wb**2*K)/LP
meanrate=sp.simplify(rate.subs(K,n*(1-q)/2))
Qblk=2*sig/LP
check(meanrate,n*Qblk,'Ensemble weighted rate = n_b Q_blk')
conditionalvar=(8*sp.pi*LP**2)**2*dt*rate/V**2
check(conditionalvar,(8*sp.pi)**2*(wl**2+r*q*wb**2*K/n)/V,
      'Conditional regional variance at fixed parities')
ensemblevar=conditionalvar.subs(K,n*(1-q)/2)
check(ensemblevar,(8*sp.pi)**2*sig/V,'Exact finite-array four-volume identity')
check(sig.subs(q,1),wl**2,'q=1 local response floor')
check(base.weights.subs(base.q,1),sp.Matrix([sp.Rational(1,2),sp.Rational(1,2),0,0]),
      'q=1 law has two equiprobable even Bell states')
rho1=base.rho.subs(base.q,1)
check(rho1,sp.diag(sp.Rational(1,2),0,0,sp.Rational(1,2)),
      'q=1 barycenter is a mixed separable state')
check(sp.trace(rho1*rho1),sp.Rational(1,2),'q=1 barycentric purity')
check(base.Q.subs(base.q,1)*base.weights.subs(base.q,1),sp.zeros(4,1),
      'q=1 full-law stationarity')
check((base.tau*sum(base.weights[k]*sum(base.state_rates[k]) for k in range(4))).subs(base.q,1),
      base.tau*base.nu/2,'q=1 mean chronogenetic activity')
# Sum exact binomial parity weights for several finite sizes; no large-N limit.
p=(1-q)/2
for nb in range(1,7):
    weighted=sum(sp.binomial(nb,k)*p**k*(1-p)**(nb-k)*conditionalvar.subs({n:nb,K:k})
                 for k in range(nb+1))
    check(weighted,((8*sp.pi)**2*sig/V).subs(n,nb),f'Finite parity mixture: {nb} block(s)')
# The discarded pathwise equality is not an identity: one generic configuration suffices.
example=(rate-meanrate).subs({n:2,K:0,q:sp.Rational(1,2),r:1,wl:1,wb:1,LP:1})
assert example!=0
passed.append('Mean relation must not be asserted for each parity realization')
report={'all_checks_passed':True,'additional_exact_checks':len(passed),
        'checks':passed,'base_exact_checks':len(base.checks),
        'scope':'Finite algebra and declared deterministic-volume, full-collection matching; no empirical data or Gaussian approximation.'}
(here/'VERIFICHE_AGGIUNTE.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
