#!/usr/bin/env python3
"""Exact auxiliary-spin inertia certificate; Python standard library only.
All physical matrices are rebuilt from spin ladder actions in rescaled Dicke
bases. No optimizer, symbolic algebra package, or Wigner symbols are used.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import json,time

def eye(d):return[[F(i==j)for j in range(d)]for i in range(d)]
def zeros(d):return[[F(0)for _ in range(d)]for _ in range(d)]
def plus(A,B,alpha=F(1)):return[[a+alpha*b for a,b in zip(x,y)]for x,y in zip(A,B)]
def scale(A,c):return[[c*x for x in row]for row in A]
def matmul(A,B):return[[sum(a*b for a,b in zip(row,col))for col in zip(*B)]for row in A]
def diagonal(xs):return[[F(x)if i==j else F(0)for j in range(len(xs))]for i,x in enumerate(xs)]

def spin_blocks(n):
 """b_r=sqrt(C(n,r))D_r, c_s=sqrt(C(3,s))D_s; total r+s conserved."""
 blocks=[]
 for total in range(n+4):
  states=[(r,total-r)for r in range(n+1)if 0<=total-r<=3]
  if not states:continue
  index={v:i for i,v in enumerate(states)};d=len(states);X=zeros(d)
  for col,(r,s)in enumerate(states):
   X[col][col]=(F(n,2)-r)*(F(3,2)-s)
   if r<n and s>0:X[index[r+1,s-1]][col]+=F((r+1)*(4-s),2)
   if r>0 and s<3:X[index[r-1,s+1]][col]+=F((n-r+1)*(s+1),2)
  G=diagonal([comb(n,r)*comb(3,s)for r,s in states]);GX=matmul(G,X)
  assert GX==list(map(list,zip(*GX)))
  blocks.append((states,X,G))
 return blocks

def spectrum(n):
 j=F(n,2);ss=F(3,2);alpha=F(n-4,2);constant=F(5*n*(n+2),16)
 total_spins=[j-ss+a for a in range(4)]
 xs=[(J*(J+1)-j*(j+1)-ss*(ss+1))/2 for J in total_spins]
 ws=[F((n+2)*(5*n-3),8),-F(n*(n-7),8),-F(3*n*n+n+6,8),-F(n*(n-7),8)]
 assert ws==[x*x-alpha*x-constant for x in xs]
 multiplicities=[F(0)]*4;projector_checks=0
 for states,X,G in spin_blocks(n):
  d=len(states);I=eye(d);W=plus(plus(matmul(X,X),X,-alpha),I,-constant)
  GW=matmul(G,W);assert GW==list(map(list,zip(*GW)))
  ps=[]
  for a,xa in enumerate(xs):
   P=I
   for b,xb in enumerate(xs):
    if a!=b:P=scale(matmul(P,plus(X,I,-xb)),1/(xa-xb))
   assert matmul(P,P)==P
   assert matmul(X,P)==scale(P,xa)
   assert matmul(W,P)==scale(P,ws[a])
   GP=matmul(G,P);assert GP==list(map(list,zip(*GP)))
   multiplicities[a]+=sum(P[i][i]for i in range(d));ps.append(P);projector_checks+=1
  summ=zeros(d)
  for P in ps:summ=plus(summ,P)
  assert summ==I
  assert all(matmul(P,Q)==zeros(d)for a,P in enumerate(ps)for b,Q in enumerate(ps)if a!=b)
 expected=[int(2*J+1)for J in total_spins];assert multiplicities==expected
 inertia=[sum(expected[a]for a,w in enumerate(ws)if w>0),sum(expected[a]for a,w in enumerate(ws)if w<0),sum(expected[a]for a,w in enumerate(ws)if w==0)]
 if n>=8:assert inertia==[n-2,3*n+6,0]
 return {'n':n,'auxiliarySpin':'3/2','alpha':str(alpha),'constant':str(constant),'totalSpins':list(map(str,total_spins)),'xEigenvalues':list(map(str,xs)),'wEigenvalues':list(map(str,ws)),'multiplicities':expected,'inertiaPositiveNegativeZero':inertia,'exactBlockProjectorChecks':projector_checks}

# Q(sqrt(11)) arithmetic for independent checks of the existing N=13 basis.
def q(a=0,b=0):return(F(a),F(b))
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def mul(x,y):return(x[0]*y[0]+11*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def smul(x,a):return(x[0]*a,x[1]*a)
def vector_add(x,y,c=F(1)):return[add(a,smul(b,c))for a,b in zip(x,y)]
def spin(v,which):
 n=len(v)-1
 if which=='z':return[smul(a,F(n,2)-r)for r,a in enumerate(v)]
 out=[q()for _ in v]
 if which=='minus':
  for r in range(n):out[r+1]=smul(v[r],r+1)
 else:
  for r in range(1,n+1):out[r-1]=smul(v[r],n-r+1)
 return out
def tensor_images(v):
 n=len(v)-1;z=spin(v,'z');lo=spin(v,'minus');hi=spin(v,'plus');c=F(n*(n+2),12)
 return[z,vector_add(spin(z,'z'),v,-c),lo,hi,vector_add(spin(lo,'z'),spin(z,'minus')),vector_add(spin(hi,'z'),spin(z,'plus')),spin(lo,'minus'),spin(hi,'plus')]
def inner(v,w,weights):
 z=q()
 for a,b,g in zip(v,w,weights):z=add(z,smul(mul(a,b),g))
 return z

def known_code_check():
 n=13;alpha=F(9,2);constant=F(975,16)
 u=[q()for _ in range(14)];u[0]=q(0,5);u[9]=q(1)
 reports=[]
 for sign in [-1,1]:
  v=[q()for _ in range(14)];v[4]=q(sign);v[13]=q(0,5);B=[u,v]
  gram=[[inner(a,b,[comb(n,r)for r in range(14)])for b in B]for a in B]
  assert gram==[[q(990),q()],[q(),q(990)]]
  for a in B:
   for b in B:
    assert all(inner(a,Tb,[comb(n,r)for r in range(14)])==q()for Tb in tensor_images(b))
  # Direct physical compression of W on C tensor auxiliary C^4.
  cols=[]
  for b in B:
   for aux in range(4):cols.append([b[r]if s==aux else q()for r in range(14)for s in range(4)])
  def Xact(x):
   out=[q()for _ in x]
   for states,X,G in spin_blocks(n):
    ids=[4*r+s for r,s in states]
    for a,idx in enumerate(ids):
     for b,jdx in enumerate(ids):out[idx]=add(out[idx],smul(x[jdx],X[a][b]))
   return out
  weights=[comb(n,r)*comb(3,s)for r in range(14)for s in range(4)]
  for x in cols:
   Xx=Xact(x);Wx=vector_add(vector_add(Xact(Xx),Xx,-alpha),x,-constant)
   assert all(inner(y,Wx,weights)==q()for y in cols)
  reports.append({'relativeEndpointPhase':'0'if sign==1 else'pi up to a basis-vector global phase','orthonormalAfterDivisionBySqrt990':True,'exactSpinTensorCrossEntries':32,'exactAncillaCompressionEntries':64})
 return reports

def verify(root=None):
 start=time.time();root=Path(root)if root else Path(__file__).resolve().parent
 # Symbolic factorization is checked as a polynomial identity: each branch is
 # quadratic in N, so four exact evaluations per branch are more than sufficient.
 for n in [8,9,13,14]:
  j=F(n,2);a=F(n-4,2);c=F(5*n*(n+2),16)
  xs=[-F(3*(n+2),4),-F(n+8,4),F(n-6,4),F(3*n,4)]
  expected=[F((n+2)*(5*n-3),8),-F(n*(n-7),8),-F(3*n*n+n+6,8),-F(n*(n-7),8)]
  assert [x*x-a*x-c for x in xs]==expected
 spectral=[spectrum(n)for n in [7,8,9,12,13,14,16,17]]
 report={'passed':True,'method':'Independent exact physical spin matrices in rescaled Dicke bases; rational spectral projectors and quadratic-field cross compressions; Python standard library only.','generalTheorem':'k_max(N,2) <= floor((N-2)/4) for every integer N>=8.','n13Conclusion':'k_max(13,2)=2.','spectralChecks':spectral,'knownN13TwoDimensionalCodes':known_code_check(),'scope':'The upper bound covers all complex symmetric subspaces, without support or symmetry assumptions. Explicit basis checks validate attainment only.','elapsedSeconds':round(time.time()-start,3)}
 (root/'verification.json').write_text(json.dumps(report,indent=2)+'\n');return report
if __name__=='__main__':
 import sys
 print(json.dumps(verify(sys.argv[1]if len(sys.argv)>1 else None),indent=2))
