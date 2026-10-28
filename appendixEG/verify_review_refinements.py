#!/usr/bin/env python3
"""Additional algebra for the targeted E/G refinements.

These checks complement (and do not replace) the proofs and prior appendix
checks. They do not establish global stochastic or nonlinear stability.
"""
from pathlib import Path
import argparse, json
import sympy as s

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--output',type=Path,default=Path(__file__).with_name('review_refinements.json'))
 args=p.parse_args();checks={}
 def zero(name,expr):
  val=s.simplify(s.expand(expr))
  if val!=0: raise AssertionError((name,val))
  checks[name]=True
 H=s.symbols('H',positive=True);tau=s.symbols('tau',real=True)
 xx=s.symbols('x1:4',real=True);q=(tau,)+xx;a=s.exp(H*tau)
 xp=a/H;xm=H*a*sum(x*x for x in xx)-1/(H*a)
 X=s.Matrix([(xp+xm)/2]+[a*x for x in xx]+[(xp-xm)/2])
 flat=s.diag(-1,1,1,1,1);g=s.diag(-1,a*a,a*a,a*a)
 dot=lambda u,v:(u.T*flat*v)[0]
 n=H*X
 zero('ambient hyperboloid',dot(X,X)-1/H**2)
 zero('spacelike unit normal',dot(n,n)-1)
 def Gamma(l,i,j):
  if l==0 and i==j and i>0:return H*a*a
  if l>0 and ((i==0 and j==l) or (j==0 and i==l)):return H
  return s.S.Zero
 for i in range(4):
  zero(f'normal tangent {i}',dot(X,X.diff(q[i])))
  for j in range(i,4):
   zero(f'induced metric {i}{j}',dot(X.diff(q[i]),X.diff(q[j]))-g[i,j])
   hess=X.diff(q[i],q[j])-sum((Gamma(l,i,j)*X.diff(q[l]) for l in range(4)),s.zeros(5,1))
   zero(f'Gauss normal coefficient {i}{j}',dot(n,hess)+H*g[i,j])
   for A in range(5):zero(f'Gauss Hessian {A}{i}{j}',hess[A]+H**2*g[i,j]*X[A])
 N,h=s.symbols('N h',positive=True);m=1-s.exp(-4*N)
 F=lambda N:8*N-3+4*s.exp(-4*N)-s.exp(-8*N)
 total=lambda N:8*N-2*(1-s.exp(-4*N))
 rem=total(N+h)-F(N)
 zero('exact remainder rearrangement',rem-(8*h+m*m-2*s.exp(-4*N)*(1-s.exp(-4*h))))
 zero('large duration remainder',s.limit(rem,N,s.oo)-(8*h+1))
 zero('time-zero boundary remainder',rem.subs(h,0)-m*m)
 S,err2,y=s.symbols('S err2 y',positive=True)
 lr=s.log(S/err2)+y*y*(1/S-1/err2)
 zero('LR derivative in variance',s.diff(lr,S)-(S-y*y)/S**2)
 # At the observed y, all physically allowed S exceed y^2.
 checks['observed residual below measurement variance']=s.Rational('0.059')**2<s.Rational('0.092')**2
 data={'passed':all(bool(v) for v in checks.values()),'checks_count':len(checks),
       'checks':{k:bool(v) for k,v in checks.items()},'sympy':s.__version__,
       'scope':'Embedding, Gauss sign, initial-data argument support, exact remainder and observed LR monotonicity.'}
 args.output.write_text(json.dumps(data,indent=2)+'\n')
 print(f"PASS: {len(checks)} refinements checks")
if __name__=='__main__':main()
