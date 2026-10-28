#!/usr/bin/env python3
"""Algebraic checks specific to the editorial refinement of WESH Section 6.

No external data and no physical matching assumptions are tested.
Run with Python 3.10+ and SymPy; output is written beside this file.
"""
from pathlib import Path
import json
import sympy as s

t=s.symbols('tau', real=True)
H=s.symbols('H', positive=True)
chi=s.Function('chi')(t)
box=-s.diff(chi,t,2)-3*H*s.diff(chi,t)
zeta=s.Rational(3,4)*(box+4*H**2*chi)
U=3*H*(H*chi-s.diff(chi,t))
V=s.diff(chi,t,2)+2*H*s.diff(chi,t)-3*H**2*chi
D=U-zeta
C00=s.diff(chi,t,2)+box/4
Csp=-H*s.diff(chi,t)-box/4
checks={}
def check(name,expr):
    checks[name]=bool(s.simplify(expr)==0)
    if not checks[name]:
        raise AssertionError(name)
check('C00_equals_U_minus_zeta',C00-D)
check('Cij_equals_D_gij_over_three',Csp-D/3)
check('homogeneous_compensating_tensor_traceless',-C00+3*Csp)
check('homogeneous_compensating_tensor_exchange',s.diff(D,t)+4*H*D+s.diff(zeta,t))
check('de_sitter_trace_filter_implies_exchange',s.diff(U,t)+4*H*U-4*H*zeta)
check('radiative_free_solution',s.diff(s.exp(-4*H*t),t)+4*H*s.exp(-4*H*t))
omega=s.symbols('omega', real=True)
check('difference_of_homogeneous_transfer_functions',4*H/(4*H-s.I*omega)-1-s.I*omega/(4*H-s.I*omega))
# Explicit MB-to-WESH relabeling in the pressureless continuity/Euler equations.
eta,k=s.symbols('eta k', real=True)
Phi=s.Function('Phi_B')(eta);Psi=s.Function('Psi_B')(eta)
phi=s.Function('phi_MB')(eta);psi=s.Function('psi_MB')(eta)
th=s.Function('theta')(eta);hh=s.Function('Hconf')(eta)
check('MB_continuity_dictionary',(-th+3*s.diff(phi,eta)).subs(phi,Psi).doit()-(-th+3*s.diff(Psi,eta)))
check('MB_Euler_dictionary',(-hh*th+k*k*psi).subs(psi,Phi)-(-hh*th+k*k*Phi))
res={'scope':'homogeneous compensating-stress identities and potential notation',
     'sympy_version':s.__version__,'groups':len(checks),'passed':sum(checks.values()),'checks':checks}
Path(__file__).with_name('section6_refinement_checks.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
