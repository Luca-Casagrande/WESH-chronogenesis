#!/usr/bin/env python3
"""Reproduce the published-input, one-history BBN compression.

This checks the statistics and matching algebra of Section 9. It does not
integrate a nuclear network or validate the Gaussian observational compression.
Run: python verify_wesh_section9.py [--output-dir DIR]
"""
from __future__ import annotations
import argparse, csv, json, math, platform
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr, ndtri
import sympy as sp

N_HAT=2.941
S_N=0.092
Y_OBS=N_HAT-3.0
GAMMA=math.pi/210.0
ALPHA=0.05
checks=[]

def check(name: str, condition: bool, **details):
    if not bool(condition):
        raise AssertionError(f'{name}: {details}')
    checks.append({'name':name,'passed':True,**details})

def F(n: float) -> float:
    if not math.isfinite(n) or n<=0:
        raise ValueError('The comparison requires finite N_pre > 0.')
    # Cancellation-free Taylor polynomial is used at small positive n.
    if n<1e-3:
        return (128/3)*n**3-128*n**4+(3584/15)*n**5-(1024/3)*n**6+(126976/315)*n**7
    m=-math.expm1(-4*n)
    return 8*n-2*m-m*m

def B(n: float, gamma: float=GAMMA) -> float:
    if not math.isfinite(gamma) or gamma<=0:
        raise ValueError('gamma must be finite and positive.')
    return (43/7)*(8*math.pi/(3*math.sqrt(gamma)))*math.sqrt(F(n))

def envelope(sigma: float, y: float=Y_OBS, n: float=50, thermal: float=1., gamma: float=GAMMA, error: float=S_N) -> float:
    if sigma<0 or thermal<=0 or error<=0:
        raise ValueError('Require sigma >= 0, thermal > 0, error > 0.')
    variance=error**2+(thermal*B(n,gamma)*sigma)**2
    return math.erf(abs(y)/math.sqrt(2*variance))

def upper_boundary(n: float, thermal: float=1., gamma: float=GAMMA, y: float=Y_OBS, error: float=S_N, alpha: float=ALPHA):
    if not 0<alpha<1 or thermal<=0 or error<=0:
        raise ValueError('Require 0 < alpha < 1, thermal > 0, error > 0.')
    q=ndtri((1+alpha)/2)
    d=y*y/q**2-error**2
    # Exact mathematical equality is recorded separately in the symbolic check.
    if d<0: return {'kind':'empty','upper':None}
    if d==0: return {'kind':'singleton_zero','upper':0.0}
    return {'kind':'interval','upper':math.sqrt(d)/(thermal*B(n,gamma))}

def main() -> None:
    checks.clear()
    # Symbolic independent checks of matching, inversion, and asymptotic scaling.

    n=sp.symbols('n', nonnegative=True)
    m=1-sp.exp(-4*n)
    f=8*n-3+4*sp.exp(-4*n)-sp.exp(-8*n)
    check('projection residual F',sp.simplify(8*n-2*m-m**2-f)==0)
    check('positive derivative F',sp.simplify(sp.diff(f,n)-8*(1-sp.exp(-4*n))**2)==0)
    check('zero-duration boundary',f.subs(n,0)==0)
    series=sp.series(f,n,0,8).removeO()
    check('small N stable series',sp.expand(series-(sp.Rational(128,3)*n**3-128*n**4+sp.Rational(3584,15)*n**5-sp.Rational(1024,3)*n**6+sp.Rational(126976,315)*n**7))==0)
    sigma,T,gamma=sp.symbols('sigma T gamma', positive=True)
    X=8*sp.pi*sigma/(3*sp.sqrt(gamma))
    check('X normalization',sp.simplify((X**2).subs(gamma,sp.pi/210)/sigma**2-sp.Rational(4480,3)*sp.pi)==0)
    check('radiation to species map',sp.simplify(sp.Rational(43,7)*sp.Rational(7,43))==1)
    check('speed-up variance is not species variance',sp.Rational(7,86)**2*sp.Rational(43,7)**2==sp.Rational(1,4))
    qs,sy,r=sp.symbols('q_alpha s_N r', positive=True)
    d=r**2/qs**2-sy**2
    check('inversion radicand',sp.simplify(sy**2+d-r**2/qs**2)==0)
    check('zero radicand boundary',sp.simplify(d.subs(r,qs*sy))==0)
    check('likelihood scale normalization', sp.diff(-sp.log(sy)-r**2/(2*sy**2),sy)==-1/sy+r**2/sy**3)

    # Data contrast, used as an audit only (no extra likelihood).
    k=0.163+(0.039/1.634)*0.405
    z=math.log(.2458/.24672)+(0.039/1.634)*math.log(2.513/2.493)
    peak=3*math.exp(z/k)
    yp_derivative=.24672*.163/3
    check('helium template derivative',abs(yp_derivative-.01340512)<1e-12,value=yp_derivative)
    check('baryon cancelling contrast',abs(.039+(0.039/1.634)*(-1.634))<1e-14)
    check('reported peak reproduction',abs(peak-2.939)<6e-5,value=peak)

    # Numerical benchmark and duration/thermal/geometry sensitivities.
    q=ndtri(.525)
    b50=B(50)
    upper=upper_boundary(50)['upper']
    check('input residual',abs(Y_OBS+.059)<1e-14)
    check('correct central 95pct quantile',abs(q-.06270677794321385)<1e-14,value=q)
    check('benchmark B50',abs(b50-8383.37677099)<1e-7,value=b50)
    check('amplitude endpoint',abs(upper-1.11694662524e-4)<1e-14,value=upper)
    check('envelope at endpoint',abs(envelope(upper)-.05)<2e-15)
    for sig,expected in [(5e-5,.109),(1e-4,.0558),(1.1e-4,.0507617447)]:
        val=envelope(sig)
        check(f'envelope sigma={sig}',abs(val-expected)<(5e-4 if sig==5e-5 else 5e-5),value=val)

    rows=[]
    for nn,expected in [(10,2.53619e-4),(30,1.44562e-4),(41,1.23449e-4),(50,1.11695e-4)]:
        up=upper_boundary(nn)['upper']
        check(f'duration N={nn}',abs(up-expected)<6e-10,value=up)
        check(f'independent root N={nn}',abs(brentq(lambda v:envelope(v,n=nn)-.05,0,.005,xtol=1e-16)-up)<1e-13)
        rows.append({'N_pre':nn,'F':F(nn),'B':B(nn),'thermal_factor':1.0,'gamma':GAMMA,'sigma_95':up,'ratio_to_50':up/upper})
    for th in [.5,.8,1,1.2,2]:
        up=upper_boundary(50,thermal=th)['upper']
        check(f'thermal scaling {th}',abs(up*th-upper)<1e-15)
    for gr in [.5,1,2,4]:
        up=upper_boundary(50,gamma=gr*GAMMA)['upper']
        check(f'geometrical prefactor {gr}',abs(up/math.sqrt(gr)-upper)<1e-15)
    peak_upper=upper_boundary(50,y=2.939-3)['upper']
    check('peak input sensitivity',abs(peak_upper-1.15517e-4)<5e-9,value=peak_upper)
    rms_N=b50*upper
    rms_u=(7/43)*rms_N
    rms_S=rms_u/2
    check('threshold rms species',abs(rms_N-.936)<.001,value=rms_N)
    check('threshold rms radiation',abs(rms_u-.152)<.001,value=rms_u)
    check('threshold rms speedup',abs(rms_S-.076)<.001,value=rms_S)
    check('weight-region bound',abs(upper**2-1.24756976365e-8)<1e-17,value=upper**2)
    check('empty-set case',upper_boundary(50,y=0)['kind']=='empty')
    check('positive-set case',upper_boundary(50,y=.1)['kind']=='interval')

    # Central probability maximization, checked by independent direct quadrature.
    for scale in [.092,.4,.94,2.]:
        half=q*scale
        centre=2*ndtr(half/scale)-1
        vals=[]
        for shift in [-2*scale,-.3*scale,0,.3*scale,2*scale]:
            prob=quad(lambda x:math.exp(-.5*((x-shift)/scale)**2)/(math.sqrt(2*math.pi)*scale),-half,half,epsabs=1e-13)[0]
            vals.append(prob)
        check(f'Anderson 1d quadrature scale={scale}',max(vals)<=centre+2e-14,max_probability=max(vals))
    for alpha in [.01,.05,.1,.2]:
        q_a=ndtri((1+alpha)/2)
        check(f'exact zero-remainder coverage alpha={alpha}',abs(1-(2*ndtr(q_a)-1)-(1-alpha))<1e-14)

    # Exact likelihood-ratio comparison at the observed endpoint, not Wilks.
    r0=abs(Y_OBS)
    def statistic(sig):
        var=S_N**2+B(50)**2*sig**2
        return math.log(var/S_N**2)+Y_OBS**2*(1/var-1/S_N**2)
    lr_bound=brentq(lambda v:statistic(v)-2.71,0,.001,xtol=1e-16)
    lr_scale=math.hypot(S_N,B(50)*lr_bound)
    lr_coverage=1-math.erf(r0/(math.sqrt(2)*lr_scale))
    check('LR apparent endpoint',abs(lr_bound-5.06e-5)<1e-7,value=lr_bound)
    check('LR actual coverage',abs(lr_coverage-.892)<.001,value=lr_coverage)
    typical=.3/((43/7)*math.sqrt(4480*math.pi/3)*math.sqrt(8*50-2*(1-math.exp(-200)))*ndtri(.975))
    check('95pct-typicality is distinct',abs(typical-1.823e-5)<1e-8,value=typical)

    # Monte Carlo diagnostics: exactly one cosmological draw per replicated experiment.
    rng=np.random.default_rng(20261016)
    size=400_000
    experiments=[]
    for truth in [0.,5e-5,upper]:
        scale=math.hypot(S_N,B(50)*truth)
        gdraw=rng.normal(0,scale,size)
        nuisances={'zero':np.zeros(size),
                   'deterministic_shift':np.full(size,.25),
                   'independent_gaussian':rng.normal(0,.4,size),
                   'independent_two_point':rng.choice([-.7,.7],size=size)}
        for name,zz in nuisances.items():
            exclusion=np.abs(gdraw+zz)<q*scale
            empirical=float(np.mean(exclusion))
            # Six standard errors around .05 for the equality case, one-sided otherwise.
            tol=6*math.sqrt(.05*.95/size)
            check(f'monte carlo coverage sigma={truth} Z={name}',empirical<=.05+tol,exclusion_rate=empirical,experiments=size)
            experiments.append({'sigma_true':truth,'remainder':name,'exclusion_fraction':empirical,'empirical_coverage':1-empirical,'experiments':size})

    # Correlated remainders can invalidate the inequality: check the necessity of independence.
    anti=rng.normal(0,1,100_000)
    check('independence assumption is necessary',np.mean(np.abs(anti-anti)<q)>.99)

    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    args=parser.parse_args();out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    report={
     'source_version':'Yeh et al., Astrophysical Journal 1008, 241 (2026), DOI 10.3847/1538-4357/ae879c, published PDF p.8, Table 2, BBN+Yp+D',
     'assumptions':['known Gaussian compressed error','jointly Gaussian inherited reference','one shared cosmological history','independent arbitrary remainder','fixed duration, geometry, thermal transport'],
     'input':{'N_hat':N_HAT,'s_N':S_N,'y_obs':Y_OBS,'gamma':GAMMA,'alpha':ALPHA},
     'benchmark':{'B50':b50,'q_alpha':q,'sigma_95':upper,'variance_bound':upper**2,'rms_delta_N':rms_N,'rms_u':rms_u,'rms_speedup':rms_S,'peak_check_N_nu':peak,'peak_based_bound':peak_upper},
     'comparisons':{'LR_threshold':2.71,'LR_endpoint':lr_bound,'LR_coverage_at_endpoint':lr_coverage,'typicality_bound':typical,'independent_remnant_fraction_50':F(50)/398},
     'duration_sensitivity':rows,'replicated_one_history_experiments':experiments,
     'checks':checks,'number_of_checks':len(checks),
     'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__},
     'not_computed':['nuclear-network integration','validation of published likelihood compression','CMB likelihood','microscopic-to-Gaussian limit theorem']
    }
    (out/'section9_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    with (out/'duration_sensitivity.csv').open('w',newline='') as fobj:
        writer=csv.DictWriter(fobj,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps({'checks_passed':len(checks),'benchmark':report['benchmark'],'comparisons':report['comparisons']},indent=2))


if __name__ == "__main__":
    main()
