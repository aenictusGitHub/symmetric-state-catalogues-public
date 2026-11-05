#!/usr/bin/env python3
"""Exact rational/6j/independent Clebsch-Gordan checks for general odd-order bounds."""
from pathlib import Path
import json
import sympy as s
from sympy.physics.wigner import wigner_6j,wigner_3j
ROOT=Path(__file__).resolve().parent

def polynomial(n,r,L,odd=True):
 return s.factor(sum(s.rf(-r,h)*s.rf(r+(3 if odd else 2),h)*s.rf(-L,h)*s.rf(L+1,h)/(s.rf(2,h)*s.rf((1 if odd else 0)-n,h)*s.rf(n+2,h)*s.factorial(h))for h in range(min(r,L)+1)))
def eigen_sixj(n,r,L,odd=True):
 j=s.Rational(n,2);jp=j-s.Rational(1,2)if odd else j
 return s.simplify(sum((n+1)*(2*q+(2 if odd else 1))*(-1)**(n+L+q)*wigner_6j(j,j,L,jp,jp,s.Rational(2*q+1,2)if odd else q)for q in range(r+1)))
def eigen_formula(n,r,L,odd=True):
 f=s.sqrt(s.Rational(n*(n+1)-L*(L+1),n*n))if odd else 1
 D=(r+1)*(r+2)if odd else(r+1)**2
 return s.simplify(D*f*polynomial(n,r,L,odd))
def eigen_direct_cg(n,r,L):
 # Independent explicit Kraus matrices: E_(ell,q)[m',m]
 # =sqrt((n+1)(2ell+1))(-1)^(j'-m') ( j' ell j; -m' q m ).
 j=s.Rational(n,2);jp=j-s.Rational(1,2)
 T=[s.sqrt(2*L+1)*(-1)**(j-m)*wigner_3j(j,L,j,-m,0,m)for m in[j-i for i in range(n+1)]]
 Tp=[s.sqrt(2*L+1)*(-1)**(jp-m)*wigner_3j(jp,L,jp,-m,0,m)for m in[jp-i for i in range(n)]]
 image=[s.S.Zero]*n;norm_checks=[]
 for a in range(r+1):
  ell=s.Rational(2*a+1,2)
  for b in range(2*a+2):
   q=ell-b;norm=s.S.Zero
   for i in range(n):
    mp=jp-i;m=mp-q
    col=j-m
    if col.is_integer and 0<=col<=n:
     v=s.sqrt((n+1)*(2*ell+1))*(-1)**(jp-mp)*wigner_3j(jp,ell,j,-mp,q,m)
     image[i]+=v*v*T[int(col)];norm+=v*v
   assert s.simplify(norm-(n+1))==0
   norm_checks.append(norm)
 return s.simplify(sum(a*b for a,b in zip(Tp,image))),len(norm_checks)
def case(n,t,k):
 odd=t%2==1;r=t//2
 assert (k*(r+1)*(r+2)==n)if odd else(k*(r+1)**2==n+1)
 roots=[L for L in range(t+1,n if odd else n+1)if polynomial(n,r,L,odd)==0]
 direct_roots=[L for L in range(t+1,n if odd else n+1)if eigen_sixj(n,r,L,odd)==0]
 assert roots==direct_roots
 no_root=not roots and(n>=7 if odd else True)
 pure_jump=bool(odd and k==1 and r>=1 and polynomial(n,r,2*r+2)!=0)
 return {'n':n,'t':t,'saturatingDimension':k,'nonforbiddenIntegerRoots':roots,'noRootSaturationExcluded':no_root,'pureOrderJumpExcluded':pure_jump,'pureEqualityPolynomialValue':str(polynomial(n,r,2*r+2))if odd and k==1 and r>=1 else None,'saturationExcluded':no_root or pure_jump,'method':'exact six-j coefficients, cross-checked by terminating rational polynomial'}
def verify():
 # The two finite coefficient identities are proved by telescoping in the note;
 # these exact checks detect indexing/normalization mistakes in the implementation.
 coeff=0
 for r in range(21):
  for h in range(r+1):
   odd=sum(2*(a+1)*s.rf(-a,h)*s.rf(a+2,h)for a in range(h,r+1))
   expected=(r+1)*(r+2)*s.rf(-r,h)*s.rf(r+3,h)/s.Rational(h+1)
   assert odd==expected
   even=sum((2*a+1)*s.rf(-a,h)*s.rf(a+1,h)for a in range(h,r+1))
   expected2=(r+1)**2*s.rf(-r,h)*s.rf(r+2,h)/s.Rational(h+1)
   assert even==expected2;coeff+=2
 checks=0
 for n in[4,5,6,7,8,11,12,13,17,18,20,24,25]:
  for r in range(min(3,(n-1)//2)+1):
   for L in range(n):
    assert s.simplify(eigen_sixj(n,r,L)-eigen_formula(n,r,L))==0,(n,r,L)
    checks+=1
 cgchecks=[]
 for n,r,orders in[(6,1,list(range(6))),(12,1,[0,1,4,10,11]),(12,2,[0,6,10]),(24,2,[0,6,10,14,18,22,23]),(20,3,[0,8,12,14,18,19])]:
  for L in orders:
   value,norms=eigen_direct_cg(n,r,L)
   assert s.simplify(value-eigen_sixj(n,r,L))==0,(n,r,L,value,eigen_sixj(n,r,L))
   cgchecks.append({'n':n,'r':r,'L':L,'exactEigenvalue':str(value),'KrausNormsChecked':norms})
 cases=[case(*a)for a in[(6,3,1),(12,3,2),(18,3,3),(24,3,4),(78,3,13),(12,5,1),(24,5,2),(20,7,1),(35,4,4)]]
 report={'passed':True,'coefficientTelescopingChecks':coeff,'allSpinSixJPolynomialComparisons':checks,'independentExplicitKrausChecks':cgchecks,'saturationCases':cases,'scope':'A root is only a necessary condition and does not establish existence. The general identities and the upper-bound proof are in odd-saturation.tex.'}
 (ROOT/'odd-saturation-verification.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({'passed':True,'sixjChecks':checks,'directKrausChecks':len(cgchecks),'cases':cases}))
 return report
if __name__=='__main__':verify()
