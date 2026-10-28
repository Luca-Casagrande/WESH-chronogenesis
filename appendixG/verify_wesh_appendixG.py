#!/usr/bin/env python3
"""Reproduce Appendix G's algebra and numerical diagnostics.

Requires Python 3, numpy, scipy, sympy. Run:
    python verify_wesh_appendixG.py --output-dir results

Checks only the declared effective calculations. No nuclear network,
microscopic CLT, or new cosmological data fit is performed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import dblquad, quad
from scipy.optimize import brentq
from scipy.special import ndtr, ndtri
import sympy as sp

CHECKS: list[dict] = []

def check(name: str, condition: bool, **details) -> None:
    if not bool(condition):
        raise AssertionError(f"{name}: {details}")
    CHECKS.append({"name": name, "passed": True, **details})

def close(a: float, b: float, atol: float = 1e-11, rtol: float = 1e-10) -> bool:
    return math.isclose(a, b, abs_tol=atol, rel_tol=rtol)

def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('results'))
    args=parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    CHECKS.clear()
    # G.1: differentiate the Friedmann residual with conserved matter.
    H,Hdot,rho,p,L,Ldot=sp.symbols('H Hdot rho p Lambda Lambda_dot', real=True)
    mu,nu=sp.symbols('mu nu', positive=True)
    E0=3*H**2-L-rho
    Ea=2*Hdot+3*H**2-L+p
    E0dot=6*H*Hdot-Ldot+3*H*(rho+p)
    check('residual derivative',sp.expand(E0dot-(3*H*(Ea-E0)-Ldot))==0)
    Hdot_mixed=sp.solve(nu*E0+mu*Ea,Hdot)[0]
    check('mixed equation residual', sp.simplify(E0dot.subs(Hdot,Hdot_mixed)+3*H*(1+nu/mu)*E0+Ldot)==0)
    check('trace damping',sp.simplify((3*(1+nu/mu)).subs(mu,3*nu))==4)
    check('radiation fraction exponent',sp.simplify((4-3*(1+nu/mu))-(1-3*nu/mu))==0)
    check('marginal exponent',sp.simplify((4-3*(1+nu/mu)).subs(mu,3*nu))==0)
    # G.2: lower-cut cone volume.
    z,zi=sp.symbols('z z_i', real=True)
    primitive=z**5/5-z**6/2+3*z**7/7-z**8/8
    cut_poly=sp.Rational(1,280)-zi**5/5+zi**6/2-3*zi**7/7+zi**8/8
    check('cone primitive',sp.expand(sp.diff(primitive,z)-z**4*(1-z)**3)==0)
    check('lower-cut polynomial',sp.simplify(sp.integrate(z**4*(1-z)**3,(z,zi,1))-cut_poly)==0)
    check('cutoff origin',cut_poly.subs(zi,0)==sp.Rational(1,280))
    check('cutoff terminal zero',cut_poly.subs(zi,1)==0)
    check('ideal cone coefficient',sp.Rational(64,3)*sp.Rational(1,280)==sp.Rational(8,105))
    for ratio in [.0001,.04,.36,.81]:
        tau=1.7;ti=tau*ratio
        # Independent proper-time quadrature; C_a cancels.
        direct=quad(lambda v: (32*math.pi/3)*v**1.5*(math.sqrt(tau)-math.sqrt(v))**3,
                    ti,tau,epsabs=1e-13,epsrel=1e-12)[0]
        analytic=64*math.pi/3*tau**4*float(cut_poly.subs(zi,math.sqrt(ratio)))
        check('cut cone proper-time quadrature',close(direct,analytic,atol=3e-12),ratio=ratio,error=abs(direct-analytic))
    # G.2: consecutive windows; x's covariance is exp(-4|v-w|) after dividing by X^2.
    n,h,v,w=sp.symbols('N h v w', positive=True)
    m=1-sp.exp(-4*n)
    future_var=32*sp.integrate(sp.integrate(sp.exp(-4*(v-w)),(w,0,v)),(v,0,h))
    check('future-window variance',sp.simplify(future_var-(8*h-2*(1-sp.exp(-4*h))))==0)
    rectangle=16*sp.integrate(sp.exp(-4*v),(v,0,n))*sp.integrate(sp.exp(-4*w),(w,0,h))
    check('past-future covariance',sp.simplify(rectangle-m*(1-sp.exp(-4*h)))==0)
    boundary_cov=4*sp.integrate(sp.exp(-4*w),(w,0,h))
    check('future orthogonality after projection',sp.simplify(rectangle-m*boundary_cov)==0)
    fn=8*n-3+4*sp.exp(-4*n)-sp.exp(-8*n)
    remainder=m**2+future_var+2*m*boundary_cov
    check('full remainder variance',sp.simplify(remainder-(8*(n+h)-2*(1-sp.exp(-4*(n+h)))-fn))==0)
    for nn,hh in [(0.2,.1),(1.,.8),(4.,2.)]:
        square=32*quad(lambda v:quad(lambda w:math.exp(-4*(v-w)),0,v,epsabs=1e-13)[0],0,hh,epsabs=1e-12)[0]
        cross=16*dblquad(lambda w,v:math.exp(-4*(w-v)),-nn,0,lambda _:0,lambda _:hh,epsabs=1e-12,epsrel=1e-12)[0]
        check('future square quadrature',close(square,8*hh-2*(1-math.exp(-4*hh))),N=nn,h=hh)
        check('cross rectangle quadrature',close(cross,(1-math.exp(-4*nn))*(1-math.exp(-4*hh))),N=nn,h=hh)
    # G.3: entropy transport, using positive symbolic variables.
    T0,TB,TS,gB,gS,s0,sB,sS=sp.symbols('T0 TB TS gB gS s0 sB sS',positive=True)
    ab_as=(TS/TB)*(sS/sB)**sp.Rational(1,3)
    energy_ratio=gB*TB**4/(gS*TS**4)
    check('thermal transfer',sp.simplify(energy_ratio*ab_as**4-(gB/gS)*(sS/sB)**sp.Rational(4,3))==0)
    aB_a0=(T0/TB)*(s0/sB)**sp.Rational(1,3)
    check('duration from entropy',sp.simplify(sp.expand_log(sp.log(aB_a0),force=True)-(sp.log(T0/TB)+sp.log(s0/sB)/3))==0)
    # G.4: published central-rate scalings, with constants as rationals.
    py,pd,qy,qd=map(sp.Rational,['.039','-1.634','.163','.405'])
    coupling=py/(-pd)
    check('baryon exponent elimination',sp.simplify(py+coupling*pd)==0)
    kN=qy+coupling*qd
    zab=math.log(.2458/.24672)+float(coupling)*math.log(2.513/2.493)
    peak=3*math.exp(zab/float(kN))
    eta10=6.120*(2.513/2.493/(peak/3)**float(qd))**(1/float(pd))
    yp_reconstructed=.24672*(eta10/6.120)**float(py)*(peak/3)**float(qy)
    d_reconstructed=2.493e-5*(eta10/6.120)**float(pd)*(peak/3)**float(qd)
    check('peak reproduction',abs(peak-2.9390322)<5e-8,value=peak)
    check('both power laws reconstructed',close(yp_reconstructed,.2458) and close(d_reconstructed,2.513e-5),eta10=eta10)
    # Statistical endpoints and actual central-tail coverage.
    gamma=math.pi/210.;N=50.;F=8*N-3+4*math.exp(-4*N)-math.exp(-8*N)
    Vtotal=8*N-2*(1-math.exp(-4*N))
    xcoeff=8*math.pi/(3*math.sqrt(gamma));b=(43/7)*xcoeff*math.sqrt(F)
    y=-.059;sn=.092;quant=ndtri(.525)
    Sigma=abs(y)/quant;sperp=math.sqrt(Sigma**2-sn**2);endpoint=sperp/b
    check('confidence boundary probability',close(2*ndtr(abs(y)/math.sqrt(sn**2+(b*endpoint)**2))-1,.05),endpoint=endpoint)
    check('endpoint digits',abs(endpoint-1.1169466252e-4)<1e-14)
    typ=.3/((43/7)*xcoeff*math.sqrt(Vtotal)*ndtri(.975))
    typ_proj=.3/(b*ndtri(.975))
    check('total typicality probability',close(2*ndtr(.3/((43/7)*xcoeff*math.sqrt(Vtotal)*typ))-1,.95))
    check('projected typicality probability',close(2*ndtr(.3/(b*typ_proj))-1,.95))
    def lr_stat(sig:float,yy:float)->float:
        ss=sn**2+(b*sig)**2
        return math.log(ss/sn**2)+yy**2*(1/ss-1/sn**2)
    lr=brentq(lambda a:lr_stat(a,y)-2.71,1e-8,.001,xtol=1e-18)
    Sd=math.sqrt(sn**2+(b*lr)**2)
    r2=(math.log(Sd**2/sn**2)-2.71)/(1/sn**2-1/Sd**2)
    check('LR rejection threshold reconstructed',close(r2,y*y))
    for yy in [0,.02,.058,.060,.3,1.,3.]:
        check('LR rejection rule', (lr_stat(lr,yy)>2.71)==(abs(yy)<abs(y)),sample=yy)
    coverage=2*(1-ndtr(abs(y)/Sd))
    rejection_quad=quad(lambda yy:math.exp(-.5*(yy/Sd)**2)/(math.sqrt(2*math.pi)*Sd),y,-y,epsabs=1e-13)[0]
    check('LR coverage quadrature',close(coverage,1-rejection_quad),coverage=coverage)
    check('LR coverage digits',abs(coverage-.8918802)<5e-8)
    # Helium Taylor diagnostic of the published constant-parameter template.
    delta=sp.symbols('Delta',real=True);yp0=sp.Rational('.24672')
    template=yp0*(1+delta/3)**qy
    linear=sp.diff(template,delta).subs(delta,0)
    quadratic=sp.diff(template,delta,2).subs(delta,0)/2
    check('helium derivative',linear==sp.Rational('.01340512'))
    check('helium second coefficient',sp.simplify(quadratic-yp0*qy*(qy-1)/18)==0)
    rmslinear=float(linear)*sperp;quadterm=float(quadratic)*sperp**2
    check('helium linear rounding',abs(rmslinear-.0125522653)<5e-11)
    check('helium quadratic rounding',abs(quadterm+.0016396371)<5e-11)
    results={
        'environment':{'python':platform.python_version(),'sympy':sp.__version__,'numpy':np.__version__,'scipy':scipy.__version__},
        'checks_count':len(CHECKS),'checks':CHECKS,
        'numbers':{'k_N':float(kN),'z_ab':zab,'peak':peak,'B50':b,'q_005':float(quant),'Sigma95':float(Sigma),'s_perp95':sperp,'sigma_endpoint':endpoint,'sigma_typ_total':float(typ),'sigma_typ_projected':float(typ_proj),'sigma_LR':lr,'coverage_LR':coverage,'helium_linear_rms':rmslinear,'helium_quadratic_at_one_sigma':quadterm},
        'scope':'Algebra and statistics of the declared effective model only; no nuclear-network validation.'}
    output=args.output_dir/'appendixG_checks.json'
    output.write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps({'passed':len(CHECKS),'output':str(output),'numbers':results['numbers']},indent=2))

if __name__=='__main__':
    main()
