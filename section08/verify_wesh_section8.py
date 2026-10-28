#!/usr/bin/env python3
"""Reproducible algebra and statistical-reference checks for WESH Section 8.

The stochastic tests sample the declared Gaussian reference limit; they do not
simulate the microscopic WESH instrument or a BBN nuclear reaction network.
Run: python verify_wesh_section8.py --output section8_checks.json
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import scipy
from scipy.integrate import quad
import sympy as sp


def run_checks() -> dict[str, Any]:
    exact: list[dict[str, Any]] = []

    def eq(name: str, actual: sp.Expr, expected: sp.Expr = sp.Integer(0)) -> None:
        residual = sp.simplify(sp.expand(actual - expected))
        if residual != 0:
            raise AssertionError(f"{name}: nonzero symbolic residual {residual}")
        exact.append({"name": name, "passed": True})

    N, h, v, w, tau, H, Ca, A, eta, yi = sp.symbols(
        "N h v w tau H Ca A eta yi", positive=True
    )
    X, sigma, Lp = sp.symbols("X sigma Lp", positive=True)
    gamma = sp.pi / 210
    m = 1 - sp.exp(-4*N)
    G = 8*N - 2*m
    F = 8*N - 3 + 4*sp.exp(-4*N) - sp.exp(-8*N)
    y = sp.symbols('y', real=True)

    beta = sp.integrate(y**4 * (1-y)**3, (y, 0, 1))
    eq("beta_5_4", beta, sp.Rational(1,280))
    volume = sp.Rational(64,3)*sp.pi*tau**4*beta
    eq("proper_time_cone_volume", volume, sp.Rational(8,105)*sp.pi*tau**4)
    eq("hubble_cone_volume", volume.subs(tau, 1/(2*H)), gamma/H**4)
    vol_eta = sp.Rational(4,3)*sp.pi*A**4*sp.integrate(v**4*(eta-v)**3,(v,0,eta))
    eq("conformal_time_cone_volume", vol_eta,
       volume.subs(tau,A*eta**2/2))
    trunc = sp.integrate(y**4*(1-y)**3, (y,yi,1))
    eq("truncated_volume_polynomial", trunc,
       sp.Rational(1,280)-yi**5/5+yi**6/2-3*yi**7/7+yi**8/8)
    eq("truncated_volume_initial_zero", trunc.subs(yi,1))
    eq("truncated_volume_ideal_limit", trunc.subs(yi,0), beta)
    eq("normalized_amplitude", (8*sp.pi)**2*sigma**2/(9*gamma),
       sp.Rational(4480,3)*sp.pi*sigma**2)

    tau1, tau2 = sp.symbols('tau1 tau2', positive=True)
    covz = (8*sp.pi)**2*sigma**2/(sp.Rational(8,105)*sp.pi*tau2**4)
    covx = covz/((3/(4*tau1**2))*(3/(4*tau2**2)))
    X2 = (8*sp.pi)**2*sigma**2/(9*gamma)
    eq("nested_cone_normalized_covariance",covx,X2*(tau1/tau2)**2)
    dn=sp.symbols('dn',nonnegative=True)
    # Direct radiation scaling tau2/tau1 = exp(2*Delta n).
    eq("radiation_scaling_to_exponential",
       ((tau1/tau2)**2).subs(tau2,tau1*sp.exp(2*dn)),sp.exp(-4*dn))

    triangle=2*sp.integrate(sp.integrate(sp.exp(-4*(v-w)),(w,0,v)),(v,0,N))
    eq("integrated_covariance_square",16*triangle,G)
    cov_ub_xb=4*sp.integrate(sp.exp(-4*(N-v)),(v,0,N))
    eq("boundary_projection_coefficient",cov_ub_xb,m)
    eq("inherited_variance_schur_complement",G-2*m*cov_ub_xb+m**2,F)
    eq("projection_future_orthogonality",
       4*sp.integrate(sp.exp(-4*(N+h-v)),(v,0,N))-m*sp.exp(-4*h))
    eq("inherited_variance_derivative",sp.diff(F,N),8*(1-sp.exp(-4*N))**2)
    eq("inherited_initial_zero",F.subs(N,0))
    eq("inherited_small_interval",sp.limit(F/N**3,N,0),sp.Rational(128,3))
    eq("inherited_large_interval_offset",sp.limit(F-8*N,N,sp.oo),-3)

    m_h=1-sp.exp(-4*h)
    eq("past_future_integrated_covariance",
       16*sp.integrate(sp.integrate(sp.exp(-4*(w-v)),(w,N,N+h)),(v,0,N)),m*m_h)
    eq("continued_total_variance",G+G.subs(N,h)+2*m*m_h,G.subs(N,N+h))
    eq("remainder_variance",G+G.subs(N,h)+2*m*m_h-F,
       G.subs(N,N+h)-F)

    # Noise of an exact OU step and of its integrated response, derived from
    # Brownian kernels rather than generated from the target covariance matrix.
    q_x=8*sp.integrate(sp.exp(-8*v),(v,0,N))
    q_xu=8*sp.integrate(sp.exp(-4*v)*(1-sp.exp(-4*v)),(v,0,N))
    q_u=8*sp.integrate((1-sp.exp(-4*v))**2,(v,0,N))
    eq("gaussian_driver_terminal_variance",q_x,1-sp.exp(-8*N))
    eq("gaussian_driver_cross_variance",q_xu,m**2)
    eq("gaussian_driver_integral_variance",q_u,F)
    eq("initial_plus_driver_integral_variance",m**2+q_u,G)
    eq("initial_plus_driver_boundary_covariance",m*sp.exp(-4*N)+q_xu,m)

    rb, rn, vb, vn, rv = sp.symbols('rb rn vb vn rv',positive=True)
    eq("nonradiative_future_factor",(1/(rv*rn*vn))/(1/(rv*rb*vb)),rb*vb/(rn*vn))
    af=sp.Function('af')(v)
    ub=sp.symbols('ub',real=True)
    ab,rhoB,rhoS,gB,gS,gsB,gsS,TB,TS = sp.symbols(
        'ab rhoB rhoS gB gS gsB gsS TB TS',positive=True)
    U=ub*rhoB*(ab/af)**4
    eq("homogeneous_radiation_transport",sp.diff(U,v)+4*sp.diff(af,v)/af*U)
    eq("thermal_common_temperature_transfer",
       (gB*TB**4/(gS*TS**4))*((TS/TB)*(gsS/gsB)**sp.Rational(1,3))**4,
       gB/gS*(gsS/gsB)**sp.Rational(4,3))
    gstar=2+sp.Rational(7,2)+3*sp.Rational(7,4)
    eq("standard_bbn_energy_count",gstar,sp.Rational(43,4))
    eq("neutrino_fraction_at_reference",sp.Rational(7,4)/gstar,sp.Rational(7,43))
    u=sp.symbols('u',real=True)
    eq("linear_speedup",sp.diff(sp.sqrt(1+u),u).subs(u,0),sp.Rational(1,2))
    eq("neutrino_speedup_at_reference",sp.Rational(7,43)/2,sp.Rational(7,86))

    def Gnum(n: float) -> float:
        return 8*n + 2*np.expm1(-4*n)

    def Fnum(n: float) -> float:
        if n < 0.05:
            return quad(lambda x:8*(-np.expm1(-4*x))**2,0,n,
                        epsabs=1e-15,epsrel=1e-12)[0]
        mm=-np.expm1(-4*n)
        return 8*n-2*mm-mm*mm

    numerics=[]
    for n in [0.001,0.01,0.1,0.5,1.,3.,10.,30.,50.]:
        integ=32*quad(lambda x:quad(lambda z:math.exp(-4*(x-z)),0,x,
                                   epsabs=1e-12,epsrel=1e-12)[0],
                      0,n,epsabs=1e-11,epsrel=1e-11)[0]
        b=4*quad(lambda x:math.exp(-4*(n-x)),0,n,
                 epsabs=1e-13,epsrel=1e-12)[0]
        f=8*quad(lambda x:(-math.expm1(-4*x))**2,0,n,
                 epsabs=1e-13,epsrel=1e-12)[0]
        if not np.isclose(integ,Gnum(n),rtol=2e-9,atol=1e-13):
            raise AssertionError(f'accumulation quadrature N={n}')
        if not np.isclose(b,-np.expm1(-4*n),rtol=2e-11,atol=1e-13):
            raise AssertionError(f'boundary quadrature N={n}')
        if not np.isclose(f,Fnum(n),rtol=2e-10,atol=1e-13):
            raise AssertionError(f'projection quadrature N={n}')
        numerics.append({'N':n,'variance_u_over_X2':float(integ),
                         'variance_inherited_over_X2':float(f),'passed':True})

    # Positive definiteness and zero covariance with a vector of future values.
    psd=[]
    for n in [0.01,0.5,2.,50.]:
        future=np.array([0.,0.13,0.7,2.1,4.])
        cov=np.empty((6,6))
        cov[0,0]=Gnum(n)
        cov[1:,1:]=np.exp(-4*np.abs(future[:,None]-future[None,:]))
        cov[0,1:]=cov[1:,0]=(-np.expm1(-4*n))*np.exp(-4*future)
        proj=np.eye(6); proj[0,1]=np.expm1(-4*n)
        cp=proj@cov@proj.T
        mineig=float(np.linalg.eigvalsh(cov).min())
        residual=float(np.max(np.abs(cp[0,1:])))
        assert mineig > -1e-11 and residual < 1e-11
        assert np.isclose(cp[0,0],Fnum(n),rtol=1e-9,atol=1e-12)
        psd.append({'N':n,'min_covariance_eigenvalue':mineig,
                    'max_future_cross_covariance':residual,'passed':True})

    rng=np.random.default_rng(16082026)
    nsamples=120000
    sampling=[]
    for n in [0.1,1.,10.,50.]:
        decay=math.exp(-4*n); mm=1-decay
        qxx=8*quad(lambda z:math.exp(-8*z),0,n)[0]
        qxu=8*quad(lambda z:math.exp(-4*z)*(1-math.exp(-4*z)),0,n)[0]
        quu=8*quad(lambda z:(1-math.exp(-4*z))**2,0,n)[0]
        L=np.linalg.cholesky(np.array([[qxx,qxu],[qxu,quu]]))
        noise=rng.standard_normal((nsamples,2))@L.T
        x0=rng.standard_normal(nsamples)
        xb=decay*x0+noise[:,0]
        ubv=mm*x0+noise[:,1]
        up=ubv-mm*xb
        var=float(np.var(up,ddof=1)); target=Fnum(n)
        zvar=(var/target-1)/math.sqrt(2/(nsamples-1))
        maxcorr=abs(float(np.corrcoef(up,xb)[0,1]))
        xf=xb.copy()
        for delta in [0.1,0.3,0.8,1.]:
            aa=math.exp(-4*delta)
            xf=aa*xf+math.sqrt(1-aa*aa)*rng.standard_normal(nsamples)
            maxcorr=max(maxcorr,abs(float(np.corrcoef(up,xf)[0,1])))
        if abs(zvar)>6 or maxcorr>6/math.sqrt(nsamples):
            raise AssertionError(f'Gaussian sample consistency N={n}')
        sampling.append({'N':n,'samples':nsamples,'variance':var,
                         'theoretical_variance':target,'variance_z_score':zvar,
                         'max_abs_correlation_with_future':maxcorr,'passed':True})

    # These factors are copied from the already approved BBN normalization;
    # they are a handoff check, not a new likelihood calculation.
    bbn=[]
    for n,rounded in [(10.,3692.),(30.,6477.),(50.,8383.)]:
        factor=43/7*math.sqrt(4480*math.pi/3*Fnum(n))
        assert abs(factor-rounded)<1.
        bbn.append({'N_pre':n,'B_of_N':factor,'matches_existing_rounded_factor':True})

    return {'scope':'Section 8, fixed-background local reference law; no nuclear-network calculation',
            'versions':{'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__},
            'symbolic_checks':exact,'symbolic_count':len(exact),
            'quadrature_checks':numerics,'positive_covariance_checks':psd,
            'gaussian_reference_sampling':sampling,'bbn_normalization_handoff':bbn,
            'all_passed':True}


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('section8_checks.json'))
    args=parser.parse_args()
    result=run_checks()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(f"PASS: {result['symbolic_count']} exact symbolic checks; "
          f"{len(result['quadrature_checks'])} quadrature cases; "
          f"{len(result['positive_covariance_checks'])} covariance cases; "
          f"{len(result['gaussian_reference_sampling'])} Gaussian-reference sampling cases.")
    print(args.output)

if __name__=='__main__':
    main()
