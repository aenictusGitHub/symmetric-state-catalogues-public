#!/usr/bin/env python3
"""Direct symmetric three-copy embedding and first-leg contraction controls."""
from pathlib import Path
from itertools import combinations_with_replacement
from collections import Counter
from math import factorial,sqrt
import json,numpy as np,sympy as s
from scipy.linalg import expm
HERE=Path(__file__).resolve().parent

def symmetric_power(c,r):
 out={}
 for w in combinations_with_replacement(range(len(c)),r):
  mult=factorial(r)
  for count in Counter(w).values():mult/=factorial(count)
  out[w]=sqrt(mult)*np.prod([c[i]for i in w])
 return out

def check(c,label,known_order=None):
 n=len(c)-1;r=3;v3=symmetric_power(c,r);v2=symmetric_power(c,r-1);pairs=list(v2)
 E=np.array([[sqrt((a.count(i)+1)/r)*v3[tuple(sorted(a+(i,)))]for i in range(n+1)]for a in pairs])
 pair=np.array(list(v2.values()));expected=pair[:,None]*c[None,:]
 error=float(np.max(abs(E-expected)));assert error<1e-14
 polys=[];checks=[]
 z=np.linspace(-1,1,n+1)
 for ell in range(7):
  p=z**ell
  for h in polys:p=p-h*np.dot(h,p)/np.dot(h,h)
  p=p/np.max(abs(p));polys.append(p)
  C=(E*p[None,:])@E.conj().T;moment=np.vdot(c,p*c)
  target=moment*np.outer(pair,pair.conj());residual=float(np.max(abs(C-target)));assert residual<1e-13
  if known_order and 1<=ell<=known_order:assert abs(moment)<1e-12
  checks.append({'rank':ell,'factorizationEntryError':residual,'normalizedDiagonalMomentAbs':float(abs(moment))})
 return {'label':label,'n':n,'copies':r,'symmetricInputDimension':len(v3),'residualPairDimension':len(v2),'embeddingMaximumEntryError':error,'contractionChecks':checks}

def main():
 cat=HERE.parents[2]/'anticoherent-state-catalogue/data/catalogue.json';data=json.loads(cat.read_text());reports=[];exactchecks=[]
 rng=np.random.default_rng(20620);c=rng.normal(size=21)+1j*rng.normal(size=21);c/=np.linalg.norm(c);reports.append(check(c,'generic complex non-AC state'))
 for stateid,order in [('acliege-n20-t5',5),('acliege-n18-t5',5),('ac-n42-t7',6)]:
  state=next(x for x in data['states']if x['id']==stateid);coeff=[s.sympify(x)for x in state['exactCoefficients']];n=state['n'];j=s.Rational(n,2);moments=[]
  for q in range(order+1):
   for p in range(order+1-q):
    value=sum(s.conjugate(coeff[i+q])*coeff[i]*(i-j)**p*s.sqrt(s.Rational(factorial(n-i)*factorial(i+q),factorial(n-i-q)*factorial(i)))for i in range(n+1-q))
    target=sum((i-j)**p for i in range(n+1))/(n+1)if q==0 else s.Integer(0)
    assert s.simplify(s.expand_complex(value-target))==0,(stateid,q,p)
    moments.append({'raisingPower':q,'diagonalPower':p,'identityVerifiedExactly':True})
  exactchecks.append({'stateId':stateid,'n':n,'orderChecked':order,'moments':moments})
  c=np.array([complex(x.evalf(25))for x in coeff]);d=n+1;z=np.diag(np.arange(d)-n/2);up=np.zeros((d,d))
  for i in range(n):up[i+1,i]=sqrt((n-i)*(i+1))
  x=(up+up.T)/2;y=(up-up.T)/2j;c=expm(-1j*(.43*x-.79*y+.31*z))@c
  reports.append(check(c,'generic rotation of exact '+stateid,order))
 report={'passed':True,'exactKnownStateMoments':exactchecks,'physicalContractionControls':reports,
 'algebraicIdentity':'The normalized symmetric embedding coefficient is sqrt((count_a(i)+1)/3). For psi^tensor3, this gives E[a,i]=psi_i*(psi^tensor2)_a. Hence contraction of the first factor with any diagonal tensor T gives <psi|T|psi>|psi^tensor2><psi^tensor2|. Rotation invariance extends this necessary condition to the twirled moment.',
 'scope':'Exact known AC5 and AC6 witness moments, plus numerical physical-embedding and contraction checks for generic complex states and generic collective rotations. These controls supplement the exact exclusion verifier; numerical errors are not used in the proof.'}
 (HERE/'positive-control.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
