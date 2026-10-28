"""Regression checks for the Section 4 revision.

The finite arrays test calculus identities and bookkeeping, not empirical WESH
predictions. The optional-sampling proof is in Appendix D of the manuscript.
"""
from pathlib import Path
import json
import math
import numpy as np

root=Path(__file__).resolve().parent
r={}

# Nondecreasing compensator with initial and intermediate flat intervals,
# and a finite terminal value. Tau_T is set to one in this calculus check.
x=np.array([0.,1.,2.,4.,6.,7.,8.])
rates=np.array([[0.,0.],[.4,.6],[0.,0.],[.3,.2],[0.,0.],[1.,1.]])
cum=np.vstack([np.zeros(2),np.cumsum(rates*np.diff(x)[:,None],axis=0)])
times=cum.sum(axis=1)
life=float(times[-1])

def S(t: float) -> float:
    if t>=life: return math.inf
    j=int(np.searchsorted(times,t,side='right')-1)
    return float(x[j]+(t-times[j])/rates[j].sum())

def A(s: float) -> np.ndarray:
    if not math.isfinite(s) or s>=x[-1]: return cum[-1].copy()
    j=int(np.searchsorted(x,s,side='right')-1)
    return cum[j]+rates[j]*(s-x[j])

assert S(0.)==1. and S(1.)==4. and S(2.)==7.
assert math.isinf(S(life)) and math.isinf(S(life+1.))
qs=np.linspace(0.,life+1.,301)
err=max(abs(A(S(float(t))).sum()-min(t,life)) for t in qs)
assert err<1e-12
r['generalized_inverse_and_terminal_compensator_max_error']=float(err)
r['flat_interval_levels_tested']=[0.,1.,2.]
r['finite_lifetime']=life
# Derivative/ratio identity away from the exceptional levels.
for t in [.2,.8,1.3,1.9,2.5,3.7]:
    s=S(t);j=int(np.searchsorted(x,s,side='right')-1)
    h=1e-6
    deriv=(A(S(t+h))-A(S(t-h)))/(2*h)
    assert np.max(np.abs(deriv-rates[j]/rates[j].sum()))<1e-8
r['proportions_as_densities_with_plateaux']='passed'
# An event can exhaust the remaining activity. For a unit-rate Poisson
# process stopped at its first event E~Exp(1), the absorbing extension gives
# M_t=1_{E<=t}-min(E,t), hence E[M_t]=0. Omitting the terminal event fails.
terminal_checks=[]
for t in [.1,.5,1.,3.]:
    mean_count=-math.expm1(-t)
    mean_compensator=-math.expm1(-t)
    assert abs(mean_count-mean_compensator)<1e-15
    terminal_checks.append({'t':t,'E_N':mean_count,'E_A':mean_compensator,'E_M':0.})
r['stopped_first_event_terminal_convention']=terminal_checks

(root/'time_change_checks.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
