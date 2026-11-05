#!/usr/bin/env python3
"""Verify the imported ACLiege formula audit using exact SymPy arithmetic.

No Julia code is evaluated. The source formulas are stored as normalized Dicke
coefficient expressions in pure-state-audit.json. Run this script with Python 3
and SymPy; NumPy is used only for the independently reconstructed sphere check.
The exact verification is separate from the geometric rendering accuracy.
"""
from pathlib import Path
from functools import reduce
from math import gcd,comb,sqrt
import json,hashlib
import sympy as s
import numpy as np
ROOT=Path(__file__).resolve().parent
DATA=ROOT/'pure-state-audit.json'

def simplify(x):return s.simplify(s.expand_complex(x))

def rho(a,b,t):
 n=len(a)-1
 return s.Matrix(t+1,t+1,lambda i,j:simplify(sum(
  a[i+l]*s.conjugate(b[j+l])*s.binomial(n-t,l)*s.sqrt(s.binomial(t,i)*s.binomial(t,j)/(s.binomial(n,i+l)*s.binomial(n,j+l)))
  for l in range(n-t+1) if a[i+l]!=0 and b[j+l]!=0)))

def vector(row):return [s.sympify(c) for c in row['exactCoefficients']]

def inner(a,b):return simplify(sum(x*s.conjugate(y) for x,y in zip(a,b)))

def partition(v):
 n=len(v)-1;support=[r for r,c in enumerate(v) if c!=0];lo,hi=min(support),max(support)
 step=reduce(gcd,(r-lo for r in support),0) or 1;y=s.Symbol('y')
 terms={int((r-lo)/step):s.simplify((-1)**r*s.sqrt(s.binomial(n,r))*v[r]) for r in support}
 base=terms[0]
 p=s.Poly(sum(s.simplify(c/base)*y**r for r,c in terms.items()),y,extension=True)
 result=[m for m in [lo,n-hi] if m]
 for factor,multiplicity in s.sqf_list(p)[1]:result.extend([int(multiplicity)]*(int(factor.degree())*step))
 return sorted(result,reverse=True)

def point_ray_error(points,v):
 n=len(v)-1;coeff=np.array([1.+0j])
 for point in points:
  x,y,z=point['x'],point['y'],point['z']
  assert abs(x*x+y*y+z*z-1)<1e-12
  a=sqrt(max(0.,(1+z)/2));b=(x+1j*y)/(2*a) if a>1e-10 else 1.+0j
  for _ in range(point['multiplicity']):coeff=np.convolve(coeff,[a,-b])
 coeff=np.array([(-1)**r*c/sqrt(comb(n,r)) for r,c in enumerate(coeff)])
 expected=np.array([complex(s.N(c,30)) for c in v]);coeff/=np.linalg.norm(coeff);expected/=np.linalg.norm(expected)
 phase=np.vdot(expected,coeff)
 return float(np.linalg.norm(expected-coeff*np.conj(phase)/abs(phase)))

def main():
 data=json.loads(DATA.read_text());pure=[];subspaces=[];ancillary=[]
 for record in data['states']:
  v=vector(record);n=record['n'];t=record['achievedOrder'];assert len(v)==n+1 and inner(v,v)==1
  for order in range(1,t+1):assert rho(v,v,order)==s.eye(order+1)/(order+1),(record['id'],order)
  next_res=rho(v,v,t+1)-s.eye(t+2)/(t+2)
  defect=simplify(sum(x*s.conjugate(x) for x in next_res));assert defect!=0
  assert s.simplify(defect-s.sympify(record['verification']['nextOrderFrobeniusDistanceSquared']))==0
  actual_partition=partition(v);assert actual_partition==record['partition']
  assert sum(p['multiplicity'] for p in record['majoranaPoints'])==n
  ray_error=point_ray_error(record['majoranaPoints'],v);assert ray_error<1e-10
  pure.append({'id':record['id'],'exactOrder':t,'partition':actual_partition,'allChecksPassed':True,'sphereReconstructionError':ray_error})
  print('Verified',record['id'],flush=True)
 for record in data['subspaces']:
  basis=list(map(vector,record['basis']));n,k,t=record['n'],record['dimension'],record['order']
  assert len(basis)==k and all(len(v)==n+1 for v in basis)
  for i,a in enumerate(basis):
   for j,b in enumerate(basis):
    assert inner(a,b)==int(i==j)
    assert rho(a,b,t)==s.eye(t+1)*int(i==j)/(t+1),(record['id'],i,j)
  subspaces.append({'id':record['id'],'dimension':k,'order':t,'orthonormal':True,'allCrossReductionsVerified':True,'scalarEntriesChecked':k*k*(t+1)**2})
  print('Verified',record['id'],flush=True)
 for record in data['ancillary']['states']:
  v=vector(record);n=record['n'];t=record['achievedOrder'];assert inner(v,v)==1 and len(v)==n+1
  if t:assert rho(v,v,t)==s.eye(t+1)/(t+1)
  assert rho(v,v,t+1)!=s.eye(t+2)/(t+2)
  sectors=[s.Integer(1)];purities={0:s.Integer(1)}
  for order in range(1,n+1):
   q=min(order,n-order)
   if q not in purities:
    r=rho(v,v,q);purities[q]=simplify(sum(x*s.conjugate(x) for x in r))
   sectors.append(simplify(s.binomial(n,order)*(2**order*purities[q]-sum(s.binomial(order,l)*sectors[l]/s.binomial(n,l) for l in range(order)))))
  assert all(s.simplify(x-s.sympify(y))==0 for x,y in zip(sectors[1:],record['sectorVectorExact']))
  ancillary.append({'id':record['id'],'exactOrder':t,'exactSectorVectorVerified':True})
 report={'auditSHA256':hashlib.sha256(DATA.read_bytes()).hexdigest(),'allChecksPassed':True,'juliaExecuted':False,
  'exactPureStates':pure,'exactSubspaces':subspaces,'ancillaryExactStates':ancillary,
  'numericalPureStatesNotPromotedToExact':[r['id'] for r in data['numericalStates']],
  'sourceCaveats':data['sourceCaveats'],'maximumSphereReconstructionError':max(r['sphereReconstructionError'] for r in pure)}
 (ROOT/'pure-state-verification.json').write_text(json.dumps(report,indent=2)+'\n')
 print('Verified',len(pure),'exact OAC formulas,',len(subspaces),'exact subspaces and',len(ancillary),'ancillary formulas.')
if __name__=='__main__':main()
