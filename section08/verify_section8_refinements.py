#!/usr/bin/env python3
"""Checks for the targeted Section 8 review and a Section 9 numerical handoff.

The observational numbers are from the published Yeh et al. ApJ 1008, 241
(2026), Table 2, BBN+Yp+D row, DOI 10.3847/1538-4357/ae879c. No nuclear network, new likelihood, or full WESH CLT is claimed.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
from scipy.stats import norm
import sympy as sp

def check_all() -> dict:
    checks=[]
    def eq(name,a,b=0):
        residual=sp.simplify(a-b)
        if residual!=0: raise AssertionError(f'{name}: {residual}')
        checks.append({'name':name,'passed':True})
    N,h=sp.symbols('N h',positive=True)
    g,sig,H1,H2,t1,t2,R,Lam,u=sp.symbols('g sig H1 H2 t1 t2 R Lam u',positive=True)
    dn=sp.symbols('dn',nonnegative=True)
    covz=(8*sp.pi)**2*sig**2/(g/H2**4)
    covx=covz/(9*H1**2*H2**2)
    eq('nested_self_similar_normalized_covariance',covx,(8*sp.pi)**2*sig**2/(9*g)*(H2/H1)**2)
    eq('radiation_exponent_independent_of_gamma',covx.subs(H2,H1*sp.exp(-2*dn)),(8*sp.pi)**2*sig**2/(9*g)*sp.exp(-4*dn))
    eq('friedmann_with_mean_lambda',(R+Lam+R*u)/(R+Lam),1+R/(R+Lam)*u)
    eq('radiation_limit_of_friedmann',((R+Lam+R*u)/(R+Lam)).subs(Lam,0),1+u)
    eq('lambda_retained_linear_speedup',sp.diff(sp.sqrt(1+R/(R+Lam)*u),u).subs(u,0),R/(2*(R+Lam)))
    m=1-sp.exp(-4*N); G=8*N-2*m
    F=8*N-3+4*sp.exp(-4*N)-sp.exp(-8*N)
    eq('variance_complement_is_boundary_projection',G-F,m**2)
    # Reduce the exactly equivalent quotient before taking the asymptotic limit.
    eq('large_N_fraction_loss',sp.limit(N*sp.factor(G-F)/G,N,sp.oo),sp.Rational(1,8))
    eq('large_N_normalized_F',sp.limit(F/N,N,sp.oo),8)
    eq('N_sensitivity_log_derivative',sp.limit(N*sp.diff(sp.log(F**(-sp.Rational(1,2))),N),N,sp.oo),-sp.Rational(1,2))
    rem=G.subs(N,N+h)-F
    eq('future_remainder_variance_positive_integral',rem,m**2+(8*h-2*(1-sp.exp(-4*h)))+2*m*(1-sp.exp(-4*h)))
    T0,TB,gs0,gsB=sp.symbols('T0 TB gs0 gsB',positive=True)
    Nth=sp.log(T0/TB)+sp.log(gs0/gsB)/3
    eq('entropy_duration_relation',sp.expand_log(sp.log((T0/TB)*(gs0/gsB)**sp.Rational(1,3)),force=True),sp.expand_log(Nth,force=True))
    bbn_x=math.sqrt(4480*math.pi/3)
    def fn(n): return 8*n-3+4*math.exp(-4*n)-math.exp(-8*n)
    bfactor=lambda n:(43/7)*bbn_x*math.sqrt(fn(n))
    N0=50.0
    measurement_mean=2.941; measurement_sd=0.092; y=measurement_mean-3
    alpha=.05
    zcentral=float(norm.ppf((1+alpha)/2))
    sdcrit=abs(y)/zcentral
    sigcrit=math.sqrt(sdcrit**2-measurement_sd**2)/bfactor(N0)
    penv=lambda s:2*norm.cdf(abs(y)/math.hypot(measurement_sd,bfactor(N0)*s))-1
    if abs(penv(sigcrit)-alpha)>1e-13: raise AssertionError('central inversion')
    if not (1.11e-4 < sigcrit < 1.13e-4): raise AssertionError('nominal benchmark changed')
    checks.append({'name':'central_tail_inversion_of_versioned_BBN_only_input','passed':True})
    checks.append({'name':'rounded_1p1e_4_benchmark_reproduced','passed':True})
    frac=fn(50)/(8*50-2*(1-math.exp(-200)))
    if abs(frac-397/398)>1e-15: raise AssertionError('retained variance fraction')
    checks.append({'name':'retained_fraction_at_50','passed':True})
    future=(8*57-2*(1-math.exp(-4*57))-fn(50))/fn(50)
    if abs(future-57/397)>1e-14: raise AssertionError('future remainder')
    checks.append({'name':'future_seven_efolds_not_boundary_remainder','passed':True})
    sensitivity=[]
    for n in [10.,30.,41.,50.]:
        sensitivity.append({'N_pre':n,'F':fn(n),'B_N':bfactor(n),
          'relative_limit_to_N50':math.sqrt(fn(50)/fn(n)),
          'sigma_95_versioned_BBN_only':math.sqrt(sdcrit**2-measurement_sd**2)/bfactor(n)})
    thermal=[]
    for T0v in [1e10,1e15]:
        for g0,gB in [(10.75,10.75),(106.75,10.75)]:
            # These counts are illustrative assignments, not a thermodynamic
            # reconstruction used in the manuscript's exact ideal model.
            n=math.log(T0v/1e-3)+math.log(g0/gB)/3
            thermal.append({'T0_GeV_assigned':T0v,'TB_GeV_assigned':1e-3,
              'gs0_assigned':g0,'gsB_assigned':gB,'N_pre':n})
    return {'scope':'Targeted algebra, second moments, and fixed-input statistical handoff only',
      'all_passed':True,'checks':checks,'check_count':len(checks),
      'boundary_fraction_N50':frac,'future_remainder_ratio_h7':future,
      'BBN_handoff':{'source':'Yeh et al., ApJ 1008, 241 (2026), Table 2, BBN+Yp+D row; DOI 10.3847/1538-4357/ae879c',
        'input_provenance':'Documented source values; the PDF is not re-read by this script','mean':measurement_mean,'measurement_sd':measurement_sd,
        'y_obs':y,'alpha':alpha,'X_per_sigma':bbn_x,'sqrt_F50':math.sqrt(fn(50)),
        'rms_u_per_sigma':bbn_x*math.sqrt(fn(50)), 'B_50':bfactor(50),
        'rms_deltaN_at_1p1e_4':bfactor(50)*1.1e-4,
        'central_probability_at_1p1e_4':float(penv(1.1e-4)),
        'sigma_95_central':sigcrit,'central_probability_at_bound':float(penv(sigcrit)),
        'sensitivity':sensitivity},'thermal_examples':thermal}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path('section8_refinements_checks.json'))
    args=p.parse_args(); result=check_all()
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(f"PASS: {result['check_count']} targeted checks")
    print(json.dumps(result['BBN_handoff'],indent=2))
