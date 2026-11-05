#!/usr/bin/env python3
from pathlib import Path
import json,hashlib
import mpmath as mp
ROOT=Path(__file__).resolve().parent
mp.mp.dps=100
source=json.loads((ROOT.parent/'analytical-constructions/t1-bases.json').read_text())
spaces={s['n']:s for s in source['subspaces']if 2<=s['n']<=20}
spaces[3]={'id':'family46-n3','n':3,'dimension':1,'order':1,'label':'GHZ line','basis':[{'exactCoefficients':['sqrt(2)/2','0','0','sqrt(2)/2']}]}
(ROOT/'subspaces.json').write_text(json.dumps(list(spaces[n]for n in range(2,21)),indent=2)+'\n')
def pair(z):return [mp.nstr(mp.re(z),85),mp.nstr(mp.im(z),85)]
def small(n):
 import sympy as sp
 basis=spaces[n]['basis'];cs=[sp.sympify(s)for s in basis[0]['exactCoefficients']];c=[mp.mpc(str(sp.N(sp.re(x),95)),str(sp.N(sp.im(x),95)))for x in cs]
 if n==2:roots=[{'root':['0','0'],'multiplicity':1},{'root':None,'multiplicity':1}]
 if n==3:roots=[{'root':pair(mp.exp(2j*mp.pi*k/3)),'multiplicity':1}for k in range(3)]
 if n==4:roots=[{'root':['0','0'],'multiplicity':1}]+[{'root':pair(mp.sqrt(2)*mp.exp(2j*mp.pi*k/3)),'multiplicity':1}for k in range(3)]
 return dict(n=n,dimension=spaces[n]['dimension'],subspaceId=spaces[n]['id'],distance='0',barycenter=['0']*3,basisCoordinates=[['1','0']]+[['0','0']]*(len(basis)-1),dickeCoefficients=list(map(pair,c)),roots=roots,multiplicityPattern=[1]*n,globalOptimalityProved=True,method='Every pure first-order AC state for N<=4 has zero Majorana barycenter; PRA90,032314(2014),SecIVB.',searchStarts=0)
results=[]
for n in range(2,21):
 if n<=4:row=small(n)
 else:
  group='5-9'if n<=9 else '10-14'if n<=14 else '15-20'
  filename='n6-radical-result.json'if n==6 else f'n{n}-result.json'
  p=ROOT/'searches'/group/filename
  if not p.exists():p=ROOT/'searches'/group/f'n{n}-candidate.json'
  row=json.loads(p.read_text())
  if n==6:
   optimized=json.loads((ROOT/'searches/5-9/n6-result.json').read_text())
   for key in ['searchStarts','searchMethods','numericalSearchMaxDistance']:row[key]=optimized[key]
   row['optimizationSource']='searches/5-9/n6-result.json'
  row['sourceFile']=str(p.relative_to(ROOT));row['globalOptimalityProved']=False
  if 'distinctRoots'in row:row['roots']=row.pop('distinctRoots')
  if n==12 and not isinstance(row['roots'][0],dict):row['roots']=[{'root':r,'multiplicity':5 if j==0 else 1}for j,r in enumerate(row['roots'])]
  if 'basisCoefficients'in row:row['basisCoordinates']=row.pop('basisCoefficients')
  if n!=6:
   root=next(r['root']for r in row['roots']if r['multiplicity']==row['dimension']);row['wronskianRoot']=root
  row['multiplicityPattern']=sorted([r['multiplicity']for r in row['roots']],reverse=True)
  row.pop('stars',None)
 row['distance']=str(row['distance']);d=mp.mpf(row['distance']);q=(mp.mpf(n-2)/n**2)**(mp.mpf(1)/(n-1));dp=abs(2*(n-1)/(n*(1+q))-1)
 if n<=4:dp=mp.mpf(0)
 row['barycentricEntanglement']=mp.nstr(1-d*d,80);row['paperDistance']=mp.nstr(dp,80);row['paperEntanglement']=mp.nstr(1-dp*dp,80)
 row['ratioToPaper']=mp.nstr(d/dp,40)if dp else None
 if not row.get('searchStarts'):row['searchStarts']=row.get('searchSummary',{}).get('totalAttempts',row.get('searchSummary',{}).get('recordedCompletedStarts',0))
 results.append(row)
 out=ROOT/'states';out.mkdir(exist_ok=True);(out/f'n{n}.json').write_text(json.dumps(row,indent=2)+'\n')
D=dict(title='Majorana barycenter displacement in selected maximal first-order AC subspaces',compiledDate='2026-11-03',scope='One specified analytical maximal first-order AC subspace at each N=2..20. Results for N>=5 are maximum candidates within those subspaces; no global maximum over all first-order AC states or all maximal subspaces is asserted.',sphereConvention='Unit stars; multiplicities counted. Root z maps to (2Re z,-2Im z,|z|²−1)/(1+|z|²); infinity maps north.',paper={'title':'Multiqubit symmetric states with maximally mixed one-qubit reductions','authors':['D.Baguette','T.Bastin','J.Martin'],'doi':'10.1103/PhysRevA.90.032314','url':'https://doi.org/10.1103/PhysRevA.90.032314','relevantSections':'SecIVB,Eqs35–37;SecIII multiplicity obstruction.'},results=results)
D['searchSummary']={'totalRecordedStartsIncludingPreviousN12':sum(r['searchStarts']for r in results),'previousN12Starts':next(r['searchStarts']for r in results if r['n']==12),'newRecordedStarts':sum(r['searchStarts']for r in results if r['n']!=12),'note':'Counts exclude discarded exploratory runs and interrupted partial batches; finite Wronskian candidates are additional.'}
(ROOT/'survey.json').write_text(json.dumps(D,indent=2)+'\n')
for r in results:print(r['n'],r['dimension'],float(r['distance']),r['multiplicityPattern'])
