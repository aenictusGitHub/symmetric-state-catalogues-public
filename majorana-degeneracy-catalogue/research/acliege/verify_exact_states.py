#!/usr/bin/env python3
"""Verify five higher-size ACLiege witnesses using exact SymPy arithmetic.

No Julia code is evaluated. The source formulas are stored as normalized Dicke
coefficient expressions in exact-states.json. Run this script with Python 3
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
DATA=ROOT/'exact-states.json'

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
    data=json.loads(DATA.read_text());pure=[]
    source=ROOT.parent.parent/data['source']
    assert hashlib.sha256(source.read_bytes()).hexdigest()==data['sourceSHA256']
    catalogue=json.loads((ROOT.parent.parent/'data/catalogue.json').read_text())
    indexed={r['id']:r for r in catalogue['states']}
    for record in data['states']:
        v=vector(record);n=record['n'];t=record['achievedOrder'];assert len(v)==n+1 and inner(v,v)==1
        for order in range(1,t+1):assert rho(v,v,order)==s.eye(order+1)/(order+1),(record['id'],order)
        next_res=rho(v,v,t+1)-s.eye(t+2)/(t+2)
        defect=simplify(sum(x*s.conjugate(x) for x in next_res));assert defect!=0
        assert s.simplify(defect-s.sympify(record['verification']['nextOrderFrobeniusDistanceSquared']))==0
        actual_partition=partition(v);assert actual_partition==record['partition']==[1]*n
        assert sum(p['multiplicity'] for p in record['majoranaPoints'])==n
        ray_error=point_ray_error(record['majoranaPoints'],v);assert ray_error<1e-10
        stored=indexed[record['id']]
        assert stored['exactCoefficients']==record['exactCoefficients']
        assert stored['partition']==record['partition'] and stored['majoranaPoints']==record['majoranaPoints']
        assert stored['verification']['exactAC'] and stored['verification']['exactMultiplicity']
        entry=next(e for e in catalogue['entries'] if e['n']==n and e['partition']==[1]*n)
        assert record['id'] in entry['stateIds'] and entry['lower']>=t
        assert entry['attainment']=='explicit' and entry['lower']<=entry['upper']
        assert entry['numericalStateIds'] and all(x.startswith('liege-') for x in entry['numericalStateIds'])
        pure.append({'id':record['id'],'exactOrder':t,'partition':actual_partition,'allChecksPassed':True,'sphereReconstructionError':ray_error,'catalogueLower':entry['lower'],'catalogueUpper':entry['upper']})
        print('Verified',record['id'],flush=True)
    assert len(pure)==5
    report={'auditSHA256':hashlib.sha256(DATA.read_bytes()).hexdigest(),'sourceSHA256':data['sourceSHA256'],'allChecksPassed':True,'juliaExecuted':False,
        'exactPureStates':pure,'maximumSphereReconstructionError':max(r['sphereReconstructionError'] for r in pure),
        'scope':'Exact AC orders and full simple-root partitions of the five attached witnesses. No upper-bound or global maximality claim follows from their next-order defects.'}
    (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Verified five exact witnesses and their catalogue lower-bound attachments.')
if __name__=='__main__':main()
