#!/usr/bin/env python3
"""Exact certificate: no pure symmetric eleven-qubit state is 4-anticoherent.
Requires Python >=3.9 and SymPy; uses only integer/rational arithmetic.
Run from any directory.  See PROOF.md for the argument.
"""
import itertools,json,math
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
HERE=Path(__file__).resolve().parent

def rank(A):return len(DomainMatrix.from_Matrix(A).convert_to(s.QQ).rref()[1])
def verify():
 data=json.loads((HERE/'certificate.json').read_text());n=11;r=3;t=4;d=12
 assert (data['n'],data['copies'],data['t'])==(n,r,t)
 states=list(itertools.combinations_with_replacement(range(d),r));weights=[2*sum(x)-r*n for x in states];occ=[tuple(x.count(i) for i in range(d)) for x in states];wi={m:[i for i,w in enumerate(weights) if w==m] for m in set(weights)}
 lower=list(itertools.combinations_with_replacement(range(d),r-1));lw=[2*sum(x)-(r-1)*n for x in lower];li={m:[i for i,w in enumerate(lw) if w==m] for m in set(lw)}
 R={};L={};indices={m:{states[idx]:k for k,idx in enumerate(src)} for m,src in wi.items()}
 for m,src in wi.items():
  if m+2 not in wi:continue
  dst=wi[m+2];A=s.zeros(len(dst),len(src));B=s.zeros(len(src),len(dst))
  for col,idx in enumerate(src):
   for k in range(n):
    if occ[idx][k]:
     ns=list(states[idx]);ns.remove(k);ns.append(k+1);row=indices[m+2][tuple(sorted(ns))];A[row,col]=occ[idx][k]*(k+1);B[col,row]=(occ[idx][k+1]+1)*(n-k)
  R[m]=A;L[m+2]=B
 blocks=[]
 for b in data['blocks']:
  J=b['J2'];H=s.Matrix([[s.Integer(x) for x in row] for row in b['highestWeight']]);mult=b['multiplicity'];assert H.shape==(len(wi[J]),mult)
  assert rank(H)==mult
  if J+2 in wi:
   assert R[J]*H==s.zeros(R[J].rows,mult)
   assert mult==len(wi[J])-rank(R[J])
  else:assert H==s.eye(len(wi[J]))
  U={J:H}
  for step,M in enumerate(range(J-2,-J-1,-2),start=1):
   U[M]=L[M+2]*U[M+2]/step
   assert all(x.q==1 for x in U[M])
   # Exact ladder relation certifies equal irrep normalization across descendants.
   assert R[M]*U[M]==(J-step+1)*U[M+2]
  den=int(b['denominatorScale']);assert den>0 and all(den%math.comb(J,k)==0 for k in range(J+1))
  blocks.append({'J2':J,'mult':mult,'U':U,'den':den})
 assert sum((b['J2']+1)*b['mult'] for b in blocks)==len(states)==364
 variables=[(ib,a,b) for ib,B in enumerate(blocks) for a in range(B['mult']) for b in range(a,B['mult'])]
 assert variables==[tuple(x) for x in data['variables']] and len(variables)==27
 def f(ell,i):
  z=2*i-n
  v=n*(n+2)
  return 1 if ell==0 else z if ell==1 else 3*z*z-v if ell==2 else z*(5*z*z-3*v+4) if ell==3 else 35*z**4-(30*v-100)*z*z+3*v*(v-8)
 rows=[]
 # Rebuild the 27 chosen necessary equations without trusting saved row entries.
 for ell in range(5):
  for other in range(ell):assert sum(f(ell,i)*f(other,i) for i in range(d))==0
 for ell,m,a,b in data['selectedConstraints']:
  assert ell in (1,2,3,4) and m in li and 0<=a<=b<len(li[m])
  pa=lower[li[m][a]];pb=lower[li[m][b]];row=[]
  for ib,alpha,beta in variables:
   B=blocks[ib];J=B['J2'];value=0
   for i in range(d):
    M=m+2*i-n
    if M not in B['U']:continue
    ua=tuple(sorted(pa+(i,)));ub=tuple(sorted(pb+(i,)));ia=indices[M][ua];ibrow=indices[M][ub];U=B['U'][M]
    product=U[ia,alpha]*U[ibrow,beta]
    if alpha!=beta:product+=U[ia,beta]*U[ibrow,alpha]
    value+=(B['den']//math.comb(J,(J-M)//2))*math.comb(n,i)*(pa.count(i)+1)*(pb.count(i)+1)*f(ell,i)*product
   row.append(int(value))
  g=math.gcd(*row);assert g>0;rows.append([x//g for x in row])
 A=s.Matrix(rows);assert A==s.Matrix([[s.Integer(x) for x in row] for row in data['selectedRows']])
 assert A.shape==(27,27) and rank(A)==27
 determinant=A.det(method='domain-ge');assert determinant!=0
 assert determinant%37==data['determinantResidue']==3 and data['determinantModulus']==37
 result={'verified':True,'statement':'No pure symmetric eleven-qubit state is 4-anticoherent.','arithmetic':'Exact integers and rational numbers only','threeCopyDimension':364,'realInvariantVariables':27,'verifiedNecessaryConstraints':27,'constraintRank':27,'kernelDimension':0,'determinantModulus':37,'determinantResidue':3}
 print(json.dumps(result,indent=2));return result
if __name__=='__main__':verify()
