#!/usr/bin/env python3
"""Exact radical association, norm certification, and all ordered cross reductions.
This verifier imports neither the constructor nor the primary verifier.  It
avoids expanding the large coefficient field: the physical reduced matrices
are computed exactly in the independently verified phase-product coordinates.
"""
from pathlib import Path
import json,hashlib
import sympy as s
HERE=Path(__file__).resolve().parent

def verify(directory=None):
 folder=Path(directory) if directory is not None else HERE
 raw=(folder/'basis.json').read_bytes();d=json.loads(raw);assert (d['n'],d['t'],d['dimension'])==(20,3,2)
 n=20;R=list(range(0,21,4));S=list(range(3,20,4));p=[s.Rational(x,24576) for x in [3927,333,9318,5018,2483,3497]];q=[s.Rational(x,384) for x in [143,23,53,117,48]]
 assert all(x>0 for x in p+q)
 for power in range(4):
  target=sum(s.Integer(r)**power for r in range(n+1))/s.Integer(n+1)
  assert sum(r**power*w for r,w in zip(R,p))==target
  assert sum(r**power*w for r,w in zip(S,q))==target
 f=[s.sqrt((n-r)*(r+1)) for r in S]
 g=[s.sqrt(s.prod((n-r-h)*(r+h+1) for h in range(3))) for r in R[:-1]]
 U=s.Rational(89947,262144);V=s.Rational(-274443,1048576);W=s.Rational(26039,393216);delta=4*U*W-V*V;assert U>0 and delta>0
 A=s.sqrt(U);B=(V+s.I*s.sqrt(delta))/(2*A)
 Z=[(-1)**i*s.binomial(4,i)*(A+i*B) for i in range(5)]
 for i,z in enumerate(Z):assert s.simplify(s.expand_complex(z*s.conjugate(z)))-f[i]**2*p[i+1]*q[i]==0
 L=[s.sqrt(s.factor(g[i]**2*p[i]*q[i])) for i in range(5)]
 # Strict rational bounds imply 234<C<237 and |L3-L4|<198<234<C<237<279<L3+L4.
 for i,(lo,hi) in enumerate([(49,50),(23,24),(261,262),(238,239),(41,42)]):assert lo*lo<L[i]**2<hi*hi
 C=-L[0]+L[1]+L[2];x=-(C*C+L[3]**2-L[4]**2)/(2*C);y=s.sqrt(L[3]**2-x*x)
 # The two remaining lengths are identities over real C, with a positive radicand by the triangle inequality just certified.
 cc,aa,bb=s.symbols('C a b',real=True);xx=-(cc*cc+aa-bb)/(2*cc);yy=aa-xx*xx
 assert s.cancel((cc+xx)**2+yy-bb)==0
 third=x+s.I*y;T=[-L[0],L[1],L[2],third,-C-third]
 # Associate the saved expressions with the complete radical recurrence exactly,
 # using canonical SymPy factor forms, rather than approximate coefficient equality.
 a=[s.sqrt(p[0])]
 for i in range(5):a.append(s.factor(g[i]*Z[i]*a[i]/(f[i]*T[i])))
 b=[s.factor(s.conjugate(Z[i]/(f[i]*a[i+1]))) for i in range(5)]
 full=[]
 for support,values,row in zip([R,S],[a,b],d['basis']):
  exact=[s.S.Zero]*(n+1)
  for r,z in zip(support,values):exact[r]=z
  saved=[s.sympify(v) for v in row['exactCoefficients']]
  for index,(actual,exported) in enumerate(zip(exact,saved)):
   assert s.factor(actual-exported)==0,('Exported radical expression differs',index,actual,exported)
  assert all(not v.free_symbols for v in saved);full.append(saved)
 # Recurrence plus |Zi|²=f_i²*p_(i+1)*q_i and |Ti|²=g_i²*p_i*q_i
 # proves |a_i|²=p_i and |b_i|²=q_i by induction from a0=sqrt(p0).
 # Supports are disjoint, so these normalized vectors are orthogonal.
 # Independently construct the physical 3-qubit reduced entries. Use formal
 # phase variables; vanishing is an identity before their radical substitution.
 alpha,beta=s.symbols('alpha beta');ell0,ell1,ell2,xx,yy=s.symbols('l0 l1 l2 x y',real=True)
 zz=[(-1)**i*s.binomial(4,i)*(alpha+i*beta) for i in range(5)]
 tt=[-ell0,ell1,ell2,xx+s.I*yy,ell0-ell1-ell2-xx-s.I*yy]
 cross={}
 for i in range(5):cross[R[i+1],S[i]]=zz[i]/f[i];cross[R[i],S[i]]=tt[i]/g[i]
 products=[{(r,r):w for r,w in zip(R,p)},cross,{(v,u):s.conjugate(z) for (u,v),z in cross.items()},{(r,r):w for r,w in zip(S,q)}]
 count=0
 for ij,P in enumerate(products):
  for aidx in range(4):
   for bidx in range(4):
    value=0
    for h in range(n-2):
     product=P.get((aidx+h,bidx+h),0)
     if product:
      value+=product*s.sqrt(s.binomial(3,aidx)*s.binomial(3,bidx))*s.binomial(n-3,h)/s.sqrt(s.binomial(n,aidx+h)*s.binomial(n,bidx+h))
    target=s.Rational(1,4) if ij in (0,3) and aidx==bidx else 0
    assert s.simplify(s.expand(value-target))==0,(ij,aidx,bidx,value);count+=1
 error=max(abs(complex(s.N(z,30))-complex(*approx)) for c,row in zip(full,d['basis']) for z,approx in zip(c,row['coefficients']));assert error<1e-12
 return {'passed':True,'n':20,'t':3,'dimension':2,'orthonormalityVerified':True,'allOrderedCrossReductionsVerified':True,'complexCrossEntriesExactlyVerified':count,'fullCrossReductionEntriesChecked':count,'exportedCoefficientRecurrenceExactlyAssociated':True,'positiveTriangleRadicandCertifiedByRationalBounds':True,'basisSHA256':hashlib.sha256(raw).hexdigest(),'maximumDisplayCoefficientError':error,'method':'Independently verify rational diagonal moments, finite-difference phase lengths and strict rational triangle bounds; associate every exported radical coefficient exactly with the norm-preserving recurrence. Then compute all 64 normalized-Dicke reduced entries symbolically in the verified phase-product coordinates. This proves all ordered reductions delta_ab I4/4 without expanding the nested radical coefficient field. No constructor or primary verifier is imported.'}
if __name__=='__main__':
 result=verify();(HERE/'independent-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
