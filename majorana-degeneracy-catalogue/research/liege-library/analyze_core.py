#!/usr/bin/env python3
"""Numerical Majorana/state checks. Observations are never exact AC proofs."""
import math
from functools import lru_cache
import mpmath as mp
mp.mp.dps=90

@lru_cache(None)
def cb(n,k):return mp.mpf(math.comb(n,k)) if 0<=k<=n else mp.mpf(0)

def normalized_points(points):
 ans=[];errors=[]
 for p in points:
  norm=mp.sqrt(sum(x*x for x in p));errors.append(abs(norm-1));ans.append(tuple(x/norm for x in p))
 return ans,max(errors,default=mp.mpf(0))

def reconstruct(points):
 """Symmetrized product of physical spinors with the supplied Bloch vectors."""
 pts,radialerror=normalized_points(points);pol=[mp.mpc(1)]
 for x,y,z in pts:
  # Choose the coordinate chart whose denominator stays away from zero.
  if z>=0:
   a=mp.sqrt((1+z)/2);b=(x+mp.j*y)/(2*a)
  else:
   b=mp.sqrt((1-z)/2);a=(x-mp.j*y)/(2*b)
  nxt=[mp.mpc(0)]*(len(pol)+1)
  for r,v in enumerate(pol):nxt[r]+=a*v;nxt[r+1]+=b*v
  pol=nxt
 n=len(points);c=[pol[r]/mp.sqrt(cb(n,r)) for r in range(n+1)];norm=mp.sqrt(sum(abs(x)**2 for x in c));c=[x/norm for x in c]
 # Fix a reproducible global phase using the largest-magnitude coefficient.
 pivot=max(range(len(c)),key=lambda i:abs(c[i]));phase=mp.conj(c[pivot])/abs(c[pivot]);c=[x*phase for x in c]
 return c,pts,radialerror

def rho_residual(c,t):
 n=len(c)-1;A=[[c[a+q]*mp.sqrt(cb(t,a)*cb(n-t,q)/cb(n,a+q)) for q in range(n-t+1)] for a in range(t+1)]
 maxerror=mp.mpf(0);frobsq=mp.mpf(0);diag=[]
 for a in range(t+1):
  for b in range(t+1):
   v=sum(A[a][q]*mp.conj(A[b][q]) for q in range(n-t+1));error=v-(mp.mpf(1)/(t+1) if a==b else 0);maxerror=max(maxerror,abs(error));frobsq+=abs(error)**2
   if a==b:diag.append(mp.re(v))
 return {'t':t,'maxEntryError':mp.nstr(maxerror,18),'frobeniusError':mp.nstr(mp.sqrt(frobsq),18),'purityExcess':mp.nstr(frobsq,18)}

def clustering(points,tolerances=(1e-12,1e-10,1e-8,1e-7,1e-6,1e-4)):
 n=len(points);dist=[[float(mp.sqrt(sum((a-b)**2 for a,b in zip(p,q)))) for q in points] for p in points];runs=[]
 for tol in tolerances:
  parent=list(range(n))
  def find(i):
   while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
   return i
  for i in range(n):
   for j in range(i):
    if dist[i][j]<=tol:parent[find(i)]=find(j)
  groups={}
  for i in range(n):groups.setdefault(find(i),[]).append(i)
  groups=sorted(groups.values(),key=lambda g:(-len(g),g[0]));within=max([dist[i][j] for g in groups for i in g for j in g]+[0.]);between=min([dist[i][j] for k,g in enumerate(groups) for h in groups[:k] for i in g for j in h]+[float('inf')])
  runs.append({'tolerance':tol,'partition':[len(g) for g in groups],'groups':groups,'maximumWithinClusterDiameter':within,'minimumInterclusterDistance':between})
 chosen=runs[3];same=[r['tolerance'] for r in runs if r['partition']==chosen['partition']]
 return {'status':'numerically inferred, not an exact degeneracy proof','selectedTolerance':chosen['tolerance'],'partition':chosen['partition'],'groups':chosen['groups'],'maxMultiplicity':max(chosen['partition']),'distinctPoints':len(chosen['partition']),'maximumWithinClusterDiameter':chosen['maximumWithinClusterDiameter'],'minimumInterclusterDistance':chosen['minimumInterclusterDistance'],'samePartitionAtTolerances':same,'sensitivity':[{'tolerance':r['tolerance'],'partition':r['partition']} for r in runs]}

def coefficient_strings(c,digits=30):return [[mp.nstr(mp.re(x),digits),mp.nstr(mp.im(x),digits)] for x in c]
