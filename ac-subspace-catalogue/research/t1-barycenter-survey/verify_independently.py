#!/usr/bin/env python3
"""Rebuild every saved state and star barycenter, independently of the searches."""
from pathlib import Path
import json,hashlib
import mpmath as mp
import sympy as sp
ROOT=Path(__file__).resolve().parent
mp.mp.dps=110

def cpx(x):return mp.mpc(*x)
def exact(x):
 z=sp.sympify(x,rational=True);return mp.mpc(str(sp.N(sp.re(z),115)),str(sp.N(sp.im(z),115)))
def point(root):
 if root is None:return [mp.mpf(0),mp.mpf(0),mp.mpf(1)]
 z=cpx(root);q=1+abs(z)**2;return [2*z.real/q,-2*z.imag/q,(abs(z)**2-1)/q]
def fmt(x):return mp.nstr(x,20)

def verify(row,source):
 n=row['n'];c=list(map(cpx,row['dickeCoefficients']));alpha=list(map(cpx,row['basisCoordinates']));k=len(alpha)
 basis=[[exact(v)for v in b['exactCoefficients']]for b in source['basis']]
 assert len(c)==n+1 and len(basis)==k
 rebuilt=[sum(alpha[a]*basis[a][r]for a in range(k))for r in range(n+1)]
 membership=max(abs(a-b)for a,b in zip(c,rebuilt));norm=abs(sum(abs(x)**2 for x in c)-1)
 A=mp.matrix(2,n)
 for a in range(2):
  for b in range(n):A[a,b]=c[a+b]*mp.sqrt(mp.binomial(n-1,b)/mp.binomial(n,a+b))
 rho=A*A.H;ac=max(abs(rho[a,b]-mp.mpf(a==b)/2)for a in range(2)for b in range(2))
 # Homogeneous polynomial factors; infinity contributes constant factor.
 polynomial=[mp.mpc(1)];points=[];count=0;infinity=0
 for rec in row['roots']:
  root=rec['root'];multiplicity=rec['multiplicity'];count+=multiplicity;points.extend([point(root)]*multiplicity)
  if root is None:infinity+=multiplicity;continue
  z=cpx(root)
  for _ in range(multiplicity):
   p=[mp.mpc(0)]*(len(polynomial)+1)
   for r,x in enumerate(polynomial):p[r]-=z*x;p[r+1]+=x
   polynomial=p
 assert count==n and len(polynomial)==n-infinity+1
 polynomial.extend([mp.mpc(0)]*infinity)
 original=[(-1)**r*mp.sqrt(mp.binomial(n,r))*c[r]for r in range(n+1)]
 pivot=max(range(n+1),key=lambda r:abs(original[r]));scale=original[pivot]/polynomial[pivot]
 reconstruction=max(abs(a*scale-b)for a,b in zip(polynomial,original))/max(abs(x)for x in original)
 bary=[sum(x[a]for x in points)/n for a in range(3)];distance=mp.sqrt(sum(x*x for x in bary));saved=mp.mpf(row['distance'])
 berror=max(abs(a-mp.mpf(b))for a,b in zip(bary,row['barycenter']));derror=abs(saved-distance)
 tolerance=mp.mpf('1e-65')
 assert max(membership,norm,ac,reconstruction,berror,derror)<tolerance,(n,list(map(fmt,[membership,norm,ac,reconstruction,berror,derror])))
 return dict(n=n,passed=True,sourceBasisMembershipError=fmt(membership),normalizationError=fmt(norm),oneQubitReductionError=fmt(ac),relativeRootProductReconstructionError=fmt(reconstruction),barycenterError=fmt(berror),distanceError=fmt(derror),recomputedDistance=mp.nstr(distance,80),starCount=count)

if __name__=='__main__':
 raw=(ROOT/'survey.json').read_bytes();data=json.loads(raw);sources={s['n']:s for s in json.loads((ROOT/'subspaces.json').read_text())}
 checks=[]
 for row in data['results']:
  checks.append(verify(row,sources[row['n']]));print(row['n'],'verified',checks[-1]['relativeRootProductReconstructionError'],flush=True)
 out=dict(passed=True,precisionDigits=110,acceptanceTolerance='1e-65',surveySHA256=hashlib.sha256(raw).hexdigest(),sourceBasisSHA256=hashlib.sha256((ROOT/'subspaces.json').read_bytes()).hexdigest(),stateCount=len(checks),checks=checks)
 (ROOT/'independent-verification.json').write_text(json.dumps(out,indent=2)+'\n')
