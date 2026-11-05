#!/usr/bin/env python3
"""Exact seven-by-nine conditional-density map for one C2 chart.

Ordinary Python+SymPy builds the polynomial matrix. Passing --kernel requires
SageMath's Python and should be externally time-limited. This is a diagnostic,
not an exclusion theorem. Odd residual and polar charts are not covered.
"""
import json,sys,time
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent
start=time.monotonic();a,b,c,z=s.symbols('a b c z',real=True)
Q=z**4+(b+s.I*c)*z*z+a
W=[[s.Poly(s.expand(Q*Q*z**(2*i)),z).nth(r)for r in range(13)]for i in range(3)]
bases=[];names=[]
for i in range(3):
 H=s.zeros(3);H[i,i]=1;bases.append(H);names.append(f'd{i}')
for i,j in [(0,1),(0,2),(1,2)]:
 H=s.zeros(3);H[i,j]=H[j,i]=1;bases.append(H);names.append(f're{i}{j}')
 H=s.zeros(3);H[i,j]=s.I;H[j,i]=-s.I;bases.append(H);names.append(f'im{i}{j}')
columns=[]
for H in bases:
 rho={}
 for aa,cc in [(0,0),(1,1),(2,2),(3,3),(0,2),(1,3)]:
  rho[aa,cc]=s.expand(sum(s.binomial(9,q)*W[i][aa+q]*s.conjugate(W[j][cc+q])*H[i,j]/(s.binomial(12,aa+q)*s.binomial(12,cc+q))for q in range(10)for i in range(3)for j in range(3)))
 columns.append([s.expand(s.binomial(3,aa)*rho[aa,aa]-rho[0,0])for aa in[1,2,3]]+[s.re(rho[0,2]),s.im(rho[0,2]),s.re(rho[1,3]),s.im(rho[1,3])])
M=s.Matrix.hstack(*[s.Matrix(v)for v in columns]);scales=[]
for row in range(7):
 scale=s.ilcm(*[s.denom(q) for value in M[row,:] for q in s.Poly(value,a,b,c).coeffs()])
 M[row,:]=M[row,:]*scale;scales.append(str(scale))
output={'scope':'Conditional-density equations for C2-even chart Q=z4+(b+ic)z2+a, a>0, R in span{1,z2,z4}.','exclusionProved':False,'variables':['a','b','c'],'densityCoordinates':names,'matrix':[[str(s.expand(v))for v in row]for row in M.tolist()],'rowScales':scales,'conditions':'M times densityCoordinates equalszero; physical residual density must be rankone positive semidefinite and nonzero.','elapsedSeconds':time.monotonic()-start}
(HERE/'analytic-c2-pencil.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({'stage':'matrix_built','shape':[7,9],'seconds':time.monotonic()-start}),flush=True)
if '--kernel' in sys.argv:
 from sage.all import PolynomialRing,QQ,matrix
 ring=PolynomialRing(QQ,names=('a','b','c'));field=ring.fraction_field()
 mat=matrix(field,[[field(str(v).replace('**','^'))for v in row]for row in M.tolist()])
 basis=mat.right_kernel().basis()
 result={'scope':output['scope'],'genericKernelDimension':len(basis),'basis':[[str(v)for v in vector]for vector in basis],'basisComplexity':[[{'numeratorDegree':int(v.numerator().total_degree()),'denominatorDegree':int(v.denominator().total_degree()),'numeratorTerms':len(v.numerator().monomials()),'denominatorTerms':len(v.denominator().monomials())}for v in vector]for vector in basis],'exactKernelVerified':all(mat*v==0 for v in basis),'elapsedSeconds':time.monotonic()-start,'exclusionProved':False}
 (HERE/'analytic-c2-kernel.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'stage':'kernel_built','dimension':len(basis),'seconds':time.monotonic()-start}),flush=True)
