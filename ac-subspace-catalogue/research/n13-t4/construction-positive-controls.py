#!/usr/bin/env python3
"""Invariant bosonic moment feasibility, with AC ideal constraints; numerical only."""
import sys,math,json,itertools,argparse,time
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize
from sympy import Rational
from sympy.physics.wigner import wigner_3j

def occupation_basis(n,r):
 a=list(itertools.combinations_with_replacement(range(n+1),r));weights=np.array([2*sum(x)-r*n for x in a]);occ=np.array([np.bincount(x,minlength=n+1) for x in a]);return a,weights,occ

def build(n,t,r):
 d=n+1; states,w,occ=occupation_basis(n,r); lower,lw,lo=occupation_basis(n,r-1)
 wi={int(m):np.flatnonzero(w==m) for m in set(w)};li={int(m):np.flatnonzero(lw==m) for m in set(lw)}
 raise_mats={}
 for m,src in wi.items():
  if m+2 not in wi:continue
  dst=wi[m+2];lookup={states[idx]:k for k,idx in enumerate(dst)};R=np.zeros((len(dst),len(src)))
  for col,idx in enumerate(src):
   for k in range(n):
    if occ[idx,k]:
     ns=list(states[idx]);ns.remove(k);ns.append(k+1);ns=tuple(sorted(ns));R[lookup[ns],col]=math.sqrt((n-k)*(k+1)*occ[idx,k]*(occ[idx,k+1]+1))
  raise_mats[m]=R
 projectors={}
 # Stabilize descent by projecting away roundoff leakage into higher-spin
 # representations; naive lowering amplifies this leakage near negative M.
 for M,idx in wi.items():
  C=raise_mats[M].T@raise_mats[M] if M in raise_mats else np.zeros((len(idx),len(idx)))
  C+=M*(M+2)/4*np.eye(len(idx))
  ev,E=np.linalg.eigh(C)
  for J in range(abs(M),r*n+1,2):
   keep=np.abs(ev-J*(J+2)/4)<1e-7
   if np.any(keep):projectors[M,J]=E[:,keep]
 blocks=[]
 for J in sorted(m for m in wi if m>=0):
  if len(wi[J])<=len(wi.get(J+2,[])):continue
  H=null_space(raise_mats[J],rcond=1e-12) if J+2 in wi else np.eye(len(wi[J]));vectors={J:H}
  for M in range(J-2,-J-1,-2):
   v=raise_mats[M].T@vectors[M+2]/math.sqrt((J+M+2)*(J-M)/4)
   E=projectors[M,J];v=E@(E.T@v)
   ev,U=np.linalg.eigh(v.T@v);v=v@(U@np.diag(1/np.sqrt(ev))@U.T)
   vectors[M]=v
  for V in vectors.values():assert np.max(abs(V.T@V-np.eye(H.shape[1])))<1e-10
  blocks.append({'J2':J,'mult':H.shape[1],'vectors':vectors})
 print('dimensions',len(states),len(lower),'blocks',[(b['J2'],b['mult']) for b in blocks],flush=True)
 ann={}
 for M,src in wi.items():
  for i in range(d):
   m=M-(2*i-n)
   if m not in li:continue
   dest=li[m];lookup={lower[idx]:k for k,idx in enumerate(dest)};A=np.zeros((len(dest),len(src)))
   for col,idx in enumerate(src):
    if occ[idx,i]:
     ns=list(states[idx]);ns.remove(i);A[lookup[tuple(ns)],col]=math.sqrt(occ[idx,i]/r)
   ann[M,i]=A
 # Hermitian tensor diagonal. Each rank's q=0 contraction suffices by covariance.
 td=[];j=Rational(n,2)
 for ell in range(1,t+1):
  td.append(np.array([math.sqrt(2*ell+1)*(-1)**(n-k)*float(wigner_3j(j,ell,j,j-k,0,k-j)) for k in range(d)]))
 # Real symmetric invariant variables. A real solution is sufficient to block proof search.
 variables=[]
 for ib,b in enumerate(blocks):
  for a in range(b['mult']):
   for c in range(a,b['mult']):variables.append((ib,a,c))
 columns=[]
 lm=sorted(li);offsets=np.cumsum([0]+[len(li[m])**2 for m in lm]); total=int(offsets[-1]);oi={m:i for i,m in enumerate(lm)}
 print('vars',len(variables),'constraints raw',t*total,flush=True)
 for nv,(ib,a,c) in enumerate(variables):
  b=blocks[ib];res=np.zeros((t,total))
  for M,V in b['vectors'].items():
   for i in range(d):
    if (M,i) not in ann:continue
    m=M-(2*i-n);z=ann[M,i]@V
    prod=np.outer(z[:,a],z[:,c]);prod=(prod+prod.T)/math.sqrt(2) if a!=c else prod
    k=oi[m];res[:,offsets[k]:offsets[k+1]]+=np.array([x[i] for x in td])[:,None]*prod.ravel()[None,:]
  columns.append(res.ravel())
  if nv%20==0:print('built',nv+1,flush=True)
 G=np.array(columns).T
 return blocks,variables,G


def main():
 import sympy as sp
 here=Path(__file__).resolve().parent
 records=[]
 for n in [12,14]:
  blocks,variables,G=build(n,4,4)
  U,s,Vh=np.linalg.svd(G,full_matrices=False)
  rank=int(np.count_nonzero(s>1e-10*s[0]));K=Vh[rank:].T
  trace=np.array([blocks[ib]['J2']+1 if a==b else 0 for ib,a,b in variables],float)
  c=np.zeros(n+1,complex)
  if n==12:c[[1,6,11]]=np.array([-np.sqrt(7),np.sqrt(11),np.sqrt(7)])/5
  else:
   source=json.loads((here.parent/'n14-t4/basis.json').read_text())
   c=np.array([complex(sp.N(sp.sympify(x),20)) for x in source['basis'][0]['exactCoefficients']])
  split=np.array([[c[a+b]*math.sqrt(math.comb(4,a)*math.comb(n-4,b)/math.comb(n,a+b)) for b in range(n-3)] for a in range(5)])
  rho_error=float(np.linalg.norm(split@split.conj().T-np.eye(5)/5));assert rho_error<1e-14
  states,w,occ=occupation_basis(n,4)
  fourth=np.array([math.sqrt(math.factorial(4)/math.prod(math.factorial(int(o)) for o in row))*np.prod(c**row) for row in occ])
  Rs=[]
  for block in blocks:
   R=np.zeros((block['mult'],block['mult']),complex)
   for M,V in block['vectors'].items():
    amp=V.T@fourth[w==M]
    R+=np.outer(amp,amp.conj())/(block['J2']+1)
   Rs.append(R)
  x=np.array([Rs[ib][a,b].real*(1 if a==b else math.sqrt(2)) for ib,a,b in variables])
  residual=float(np.linalg.norm(G@x));tr=float(trace@x)
  print('CONTROL',n,'residual',residual,'trace',tr,'null',K.shape[1],'singular',s[-5:],flush=True)
  assert residual<1e-8 and abs(tr-1)<1e-8 and K.shape[1]>0
  assert np.linalg.norm(K.T@trace)>1e-6
  assert max(np.linalg.norm(R.real-R.real.T) for R in Rs)<1e-12
  assert min(np.linalg.eigvalsh(R.real).min() for R in Rs)>-1e-12
  records.append({'n':n,'t':4,'copies':4,'variables':len(variables),'rank':rank,'nullDimension':int(K.shape[1]),'knownExactWitnessTrace':tr,'directFourthReductionDeviation':rho_error,'knownExactWitnessConstraintNorm':residual,'realTwirlBlocksPositiveSemidefinite':True,'smallestSingularValues':s[-10:].tolist(),'traceProjectionToNullNorm':float(np.linalg.norm(K.T@trace))})
 result={'passed':True,'stabilization':'Each lowering step is projected to the proper total-spin Casimir eigenspace before polar reorthonormalization, preventing amplified higher-spin roundoff leakage.', 'rankThresholdRelative':1e-10,'scope':'Numerical positive controls of the identical physical four-copy contraction construction. Exact known witnesses survive and have nonzero invariant trace; this verifies the scheme does not erroneously exclude adjacent known cases. This is not an exact exclusion proof.','records':records}
 (here/'construction-positive-controls.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
