#!/usr/bin/env python3
"""Exact rational certificate of Majorana partition (2,1^9) on the root box.
Uses only the Python standard library. No numerical root clustering.
"""
import json,math
from pathlib import Path
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent
Z=(Q(0),Q(0));ONE=(Q(1),Q(0))
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def scale(a,r):return (a[0]*r,a[1]*r)
def size(a):return abs(a[0])+abs(a[1])
def verify():
 d=json.loads((HERE/'witness.json').read_text());cert=json.loads((HERE/'majorana-certificate.json').read_text());center=list(map(Q,d['allCoordinateCenter']));radius=Q(d['radius']);a=[scale((center[r],center[r+10]),(-1)**r*math.comb(11,r)) for r in range(10)];b=[scale(a[r+1],r+1) for r in range(9)];S=[[Z]*17 for _ in range(17)]
 for i in range(8):
  for j,c in enumerate(reversed(a)):S[i][i+j]=c
 for i in range(9):
  for j,c in enumerate(reversed(b)):S[i+8][i+j]=c
 R=[[tuple(map(Q,z)) for z in row] for row in cert['inverseApproximation']];assert len(R)==17 and all(len(row)==17 for row in R)
 residual=Q(0)
 for i in range(17):
  rownorm=Q(0)
  for j in range(17):
   z=ONE if i==j else Z
   for k in range(17):z=add(z,scale(mul(R[i][k],S[k][j]),-1))
   rownorm+=size(z)
  residual=max(residual,rownorm)
 invnorm=max(sum(size(z) for z in row) for row in R)
 # Both real and imaginary coordinate errors are at most radius; this is
 # conservative for fixed coordinates. Derivative rows dominate polynomial rows.
 variation=2*radius*sum(r*math.comb(11,r) for r in range(10))
 bound=residual+invnorm*variation
 assert bound<1 and abs(center[9])>radius
 assert cert['partition']==[2]+[1]*9
 result={'verified':True,'method':'Exact rational bound ||I-R S||∞<1 for every polynomial in the algebraic-state box','partition':[2]+[1]*9,'maxMultiplicity':2,'distinctPoints':10,'finiteDegree':9,'infinityMultiplicity':2,'finiteRootsAllSimple':True,'inversePerturbationBoundApprox':float(bound)}
 print(json.dumps(result,indent=2));return result
if __name__=='__main__':verify()
