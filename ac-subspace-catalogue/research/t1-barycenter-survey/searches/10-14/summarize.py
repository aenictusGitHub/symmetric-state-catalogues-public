from survey import *
summary=[];base=[]
for n in [10,11,13,14]:
 s,V,M,P=get(n);k=s['dimension'];candidate=json.loads((ROOT/f'n{n}-result.json').read_text());w=json.loads((ROOT/f'n{n}-wronskian.json').read_text());search=json.loads((ROOT/f'n{n}-search.json').read_text());poles=json.loads((ROOT/f'n{n}-poles.json').read_text());lower=json.loads((ROOT/f'n{n}-lowerfold.json').read_text())
 for row in poles['attempts']:
  support=[r for r in range(n+1)if any(P[r][j]!=0 for j in row['active'])];row['zeroMultiplicity']=min(support);row['infinityMultiplicity']=n-max(support)
 (ROOT/f'n{n}-poles.json').write_text(json.dumps(poles,indent=2))
 for j in range(k):
  support=[r for r in range(n+1)if P[r][j]!=0];lo,hi=min(support),max(support);step=math.gcd(*(r-lo for r in support)) or 1
  q=[P[r][j]for r in range(lo,hi+1,step)];ys=mp.polyroots(q[::-1],maxsteps=3000,extraprec=100) if len(q)>1 else []
  roots=[mp.mpc(0)]*lo
  for y in ys:
   roots.extend([mp.root(y,step,a)for a in range(step)])
  xyz,b,d=points_bary(roots,n,n-hi);base.append(dict(n=n,basisIndex=j,distance=mp.nstr(d,70),barycenter=list(map(lambda x:mp.nstr(x,70),b))))
 candidate['searchSummary']={'unrestrictedBFGS':220,'unrestrictedDECharts':k,'fixedPoleSearchesAndBasisEndpoints':len(poles['attempts']),'lowerMultiplicityManifoldDE':8,'finiteWronskianCandidates':len(w['candidates']),'totalAttempts':len(search['attempts'])+len(poles['attempts'])+8,'unrestrictedBest':search['best']['distance'],'fixedPoleBest':poles['best']['distance'],'lowerMultiplicityNote':'Double-precision apparent excesses were root-conditioning artifacts. Independent100-digit evaluations of all eight lowerfold final roots were below the high-precision Wronskian candidate.'}
 # Fixed south pole at N14 follows identically from support, not numerical clustering.
 fixedzero=min(r for r in range(n+1)if any(P[r][j]!=0 for j in range(k)))
 if fixedzero:
  simple=min((x for x in candidate['distinctRoots']if x['root']is not None and x['multiplicity']==1),key=lambda x:abs(mp.mpc(*x['root'])))
  assert abs(mp.mpc(*simple['root']))<mp.mpf('1e-75');simple['root']=['0','0'];simple['evidence']='Exact common monomial factor of all basis polynomials.'
 candidate['sourceBasis']=[v['exactCoefficients']for v in s['basis']]
 (ROOT/f'n{n}-result.json').write_text(json.dumps(candidate,indent=2))
 summary.append({key:candidate[key]for key in ['n','dimension','subspaceId','distance','barycenter','searchSummary','optimalityProved']})
(ROOT/'basis-baseline.json').write_text(json.dumps(base,indent=2));(ROOT/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
