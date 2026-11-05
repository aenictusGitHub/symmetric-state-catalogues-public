#!/usr/bin/env python3
"""Independent exact Dicke-reduction and rational existence verification.
No primary verifier, constructor, optimizer or source symmetry ansatz is used.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import json,hashlib
import sympy as s
HERE=Path(__file__).resolve().parent

def verify():
 c=json.loads((HERE/'certificate.json').read_text());b=json.loads((HERE/'basis.json').read_text())
 n,t,k=b['n'],b['t'],b['dimension'];assert(n,t,k)==(24,2,5)
 unknown=c['unknownSlotIndices'];fixed=c['fixedSlotIndices'];m=len(unknown);total=len(c['slots'])
 assert m==23 and total==26 and len(set(unknown+fixed))==26 and set(unknown+fixed)==set(range(26))
 D=int(c['centerDenominator']);Q=int(c['inverseDenominator']);r=F(c['radius'])
 a=list(map(int,c['fullCenterNumerators']));A=[list(map(int,row)) for row in c['inverseNumerators']]
 assert len(a)==total and len(A)==m and all(len(row)==m for row in A) and D>0 and Q>0 and r>0
 fixedvals={i:F(int(c['fixedNumerators'][str(i)]),int(c['fixedDenominator'])) for i in fixed}
 assert all(F(a[i],D)==v for i,v in fixedvals.items())
 names=b['algebraicDefinition']['variables'];xs=s.symbols(' '.join(names),real=True);local=dict(zip(names,xs))
 atoms={slot:xs[j] for j,slot in enumerate(unknown)};atoms.update({i:s.Rational(v) for i,v in fixedvals.items()})
 intervals=b['algebraicDefinition']['parameterIntervals'];assert set(intervals)==set(names)
 center={xs[j]:s.Rational(a[i],D) for j,i in enumerate(unknown)}
 for j,i in enumerate(unknown):assert list(map(F,intervals[names[j]]))==[F(a[i],D)-r,F(a[i],D)+r]
 equations=c['equations'];assert len(equations)==m
 fs=[];fv=[];jac=[]
 for e in equations:
  poly=s.Integer(e['constant']);value=e['constant']*D*D;grad=[0]*total
  for i,j,w in e['terms']:
   poly+=w*atoms[i]*atoms[j];value+=w*a[i]*a[j];grad[i]+=w*a[j];grad[j]+=w*a[i]
  fs.append(s.expand(poly));fv.append(value);jac.append([grad[i] for i in unknown])
 assert all(s.expand(s.sympify(expr,locals=local)-f)==0 for expr,f in zip(b['algebraicDefinition']['equations'],fs))
 # Verify the physical reduced matrices directly, without using moment equations
 # or the constructor's coefficient map. c_r=sqrt(C(n,r))*d_r cancels the
 # denominator in Dicke splitting; only rational polynomial coefficients remain.
 amps=[];maximum_error=0
 for row in b['basis']:
  cs=[s.sympify(z,locals=local) for z in row['exactCoefficients']]
  assert len(cs)==n+1 and all(z.free_symbols<=set(xs) for z in cs)
  amps.append([s.simplify(z/s.sqrt(comb(n,i))) for i,z in enumerate(cs)])
  for z,saved in zip(cs,row['coefficients']):
   error=abs(complex(s.N(z.subs(center),30))-complex(*saved));maximum_error=max(maximum_error,error);assert error<3e-12
 physical=[]
 for i in range(k):
  for j in range(k):
   gram=sum(comb(n,r)*amps[i][r]*s.conjugate(amps[j][r]) for r in range(n+1))-int(i==j)
   physical.extend([s.expand(s.re(gram)),s.expand(s.im(gram))])
   for aa in range(t+1):
    for bb in range(t+1):
     poly=sum(comb(n-t,d)*amps[i][aa+d]*s.conjugate(amps[j][bb+d]) for d in range(n-t+1))
     if i==j and aa==bb:poly-=s.Rational(1,(t+1)*comb(t,aa))
     physical.extend([s.expand(s.re(poly)),s.expand(s.im(poly))])
 # Rational row-space inclusion proves each full reduction identity at every
 # zero of F, rather than only at the numerical center or on selected moments.
 polys=[s.Poly(p,*xs,domain=s.QQ) for p in fs+physical]
 monomials=sorted(set().union(*(set(p.monoms()) for p in polys)))
 rows=[[p.coeff_monomial(mon) for mon in monomials] for p in polys]
 R,pivots=s.Matrix(rows[:m]).rref();assert len(pivots)==m
 for row in rows[m:]:
  assert all(row[col]==sum(row[pivots[j]]*R[j,col] for j in range(m)) for col in range(len(monomials)))
 # Banach's fixed-point theorem for x -> x-A F(x), checked with integers.
 defect=F(max(sum(abs((Q*D if i==j else 0)-sum(A[i][p]*jac[p][j] for p in range(m))) for j in range(m)) for i in range(m)),Q*D)
 correction=F(max(abs(sum(A[i][p]*fv[p] for p in range(m))) for i in range(m)),Q*D*D)
 inverse_norm=F(max(sum(abs(v) for v in row) for row in A),Q)
 active=set(unknown);h=max(2*sum(abs(w) for i,j,w in e['terms'] if i in active and j in active) for e in equations)
 contraction=defect+inverse_norm*h*r;selfmap=correction+contraction*r
 assert defect<1 and contraction<1 and selfmap<r
 saved=c['bounds']
 assert F(saved['inverseDefect'])==defect and F(saved['newtonCorrection'])==correction
 assert F(saved['inverseInfinityNorm'])==inverse_norm and int(saved['jacobianLipschitzBound'])==h
 assert F(saved['contractionFactor'])==contraction and F(saved['selfMapRadiusBound'])==selfmap
 return {'passed':True,'n':n,'t':t,'dimension':k,'maximumDimension':5,'realVariables':m,'coefficientExpressionsAssociated':k*(n+1),'orthonormalityEntriesChecked':k*k,'complexCrossEntriesExactlyVerified':k*k*(t+1)**2,'fullCrossReductionEntriesChecked':225,'physicalPolynomialRowSpanVerified':True,'strictContraction':True,'strictSelfMap':True,'radius':str(r),'contractionFactor':str(contraction),'selfMapBound':str(selfmap),'certificateSHA256':hashlib.sha256((HERE/'certificate.json').read_bytes()).hexdigest(),'basisSHA256':hashlib.sha256((HERE/'basis.json').read_bytes()).hexdigest(),'maximumDisplayCoefficientError':maximum_error,'method':'Reconstruct all ordered Dicke cross reductions directly from the exported coefficient formulas, prove exact rational row-space inclusion in the root equations, independently recompute all rational contraction inequalities and isolation intervals. No builder or primary verifier imported.','scope':'Exact isolated algebraic-root basis; no short radical formula or exhaustive classification claimed.'}

if __name__=='__main__':
 report=verify();(HERE/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
