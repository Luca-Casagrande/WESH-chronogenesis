#!/usr/bin/env python3
"""Independent refinement checks for published-source matching and helium curvature.

Run next to verify_wesh_section9.py. Optionally pass the supplied published
PDF with --published-pdf PATH to reproduce the provenance checks (PyMuPDF).
This tests source extraction, algebra and numerical diagnostics, not a nuclear
reaction network or the validity of the Gaussian sampling approximation.
"""
from __future__ import annotations
import argparse, hashlib, json, math, re
from pathlib import Path
import sympy as sp
from scipy.special import ndtr
import verify_wesh_section9 as base


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--published-pdf',type=Path)
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    checks=[]
    def check(name,condition,**values):
        if not bool(condition): raise AssertionError((name,values))
        checks.append(dict(name=name,passed=True,**values))

    # Source formula algebra, at the published central neutron lifetime and rates.
    a=sp.Rational(39,1000);b=-sp.Rational(1634,1000)
    c=sp.Rational(163,1000);d=sp.Rational(405,1000)
    eta,N=sp.symbols('eta N',positive=True)
    eta0=sp.Rational(6120,1000); Y0=sp.Rational(24672,100000)
    D0=sp.Rational(2493,1000)*sp.Rational(1,100000)
    Y=Y0*(eta/eta0)**a*(N/3)**c
    DH=D0*(eta/eta0)**b*(N/3)**d
    check('baryon exponent cancels exactly',a+(a/(-b))*b==0)
    check('contrast neutrino coefficient',c+(a/(-b))*d==sp.Rational(282137,1634000))
    z=math.log(.2458/float(Y0))+(float(a)/float(-b))*math.log(2.513e-5/float(D0))
    k=float(c+(a/(-b))*d)
    check('evaluated contrast',abs(z-(-.0035451786121970747))<2e-16,value=z)
    peak=3*math.exp(z/k)
    check('evaluated peak from both published abundances',abs(peak-2.939032195479)<2e-12,value=peak)
    eta_fit=float(eta0)*((2.513e-5/float(D0))/(peak/3)**float(d))**(1/float(b))
    check('helium forward substitution',abs(float(Y.subs({eta:eta_fit,N:peak}))-.2458)<1e-13)
    check('deuterium forward substitution',abs(float(DH.subs({eta:eta_fit,N:peak}))-2.513e-5)<1e-17)

    YN=Y.subs(eta,eta0)
    d1=sp.diff(YN,N).subs(N,3);half_d2=sp.diff(YN,N,2).subs(N,3)/2
    check('exact helium first derivative',d1==sp.Rational(1340512,100000000))
    check('exact helium quadratic coefficient',half_d2==Y0*c*(c-1)/18)
    check('helium curvature negative',half_d2<0)
    sig=base.upper_boundary(50)['upper'];width=base.B(50)*sig
    rms_linear=float(d1)*width
    quadratic=float(half_d2)*width**2
    ratio=rms_linear/.0013
    check('endpoint retains published BBN only inputs',base.N_HAT==2.941 and base.S_N==.092)
    check('source width unchanged',abs(width-.9363784392437169)<1e-13,value=width)
    check('linear helium rms from inherited width',abs(rms_linear-.0125522653434747)<1e-13,value=rms_linear)
    check('linear helium relative to measurement width',abs(ratio-9.6555887257498)<1e-10,value=ratio)
    check('quadratic helium term at one inherited width',abs(quadratic+.00163963705306576)<1e-13,value=quadratic)
    exact_residual={}
    for sign in [-1,1]:
        delta=sign*width
        pred=float(YN.subs(N,3+delta))
        linear=float(Y0)+float(d1)*delta
        res=pred-linear; exact_residual[str(sign)]=res
        check(f'nonlinear constant-template residual sign={sign}',res<0,value=res)
    # These checks do not extrapolate the nuclear template over Gaussian tails.
    check('one-sigma template points physically inside Nnu>0',3-width>0)
    check('endpoint central probability unchanged',abs(base.envelope(sig)-.05)<1e-14)
    check('confidence and typicality differ quantitatively',sig/1.8235089369281727e-5>6)
    # Published restricted extra-species bounds are not dispersion bounds.
    check('BBN+CMB input not used',base.N_HAT!=2.925 and base.S_N!=.082)
    source=None
    if args.published_pdf:
        import fitz
        doc=fitz.open(args.published_pdf)
        check('published PDF has 12 pages',len(doc)==12)
        page1=doc[0].get_text(); p3=doc[2].get_text();p8=doc[7].get_text()
        check('publisher DOI in header','10.3847/1538-4357/ae879c' in page1)
        check('published table includes mean and width',all(x in p8 for x in ['2.941','0.092','2.939']))
        check('published table removes editorial query','to be checked' not in p8.lower())
        check('published response scalings match constants',all(x in p3 for x in ['0.24672','2.493','0.163','0.405','1.634']))
        source={'filename':args.published_pdf.name,'sha256':hashlib.sha256(args.published_pdf.read_bytes()).hexdigest(),'pages':12,'input_table_page':8,'scaling_equations_page':3}
    result={'checks_passed':len(checks),'checks':checks,'source':source,
            'numerics':{'z_ab':z,'k_N':k,'N_nu_central_rate':peak,'eta10_central_rate':eta_fit,
             'helium_derivative':float(d1),'helium_half_second_derivative':float(half_d2),
             'sigma_upper_95':sig,'inherited_width':width,'helium_linear_rms':rms_linear,
             'helium_linear_in_observational_sigmas':ratio,'helium_quadratic_at_one_width':quadratic,
             'helium_exact_minus_linear_at_one_width':exact_residual},
            'limitations':['constant-parameter helium power law is a diagnostic, not a new network calculation',
             'Gaussian compression and history-template validity remain model assumptions']}
    (args.output_dir/'section9_refinement_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'checks_passed':len(checks),'numerics':result['numerics']},indent=2))

if __name__=='__main__': main()
