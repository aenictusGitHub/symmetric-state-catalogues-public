from pathlib import Path
import json,re
ROOT=Path(__file__).resolve().parent
read=lambda p:json.loads((ROOT/p).read_text())
prior_poles={n:[]for n in range(15,21)}
for line in (ROOT/'poles.log').read_text().splitlines():
 match=re.fullmatch(r'(\d+) (zero|infinity) (\d+) (\d+) ([0-9.eE+-]+)',line)
 if match:
  n=int(match[1]);l=int(match[4]);prior_poles[n].append({'pole':match[2],'constraintMultiplicity':int(match[3]),'nullDimension':l,'distance':float(match[5]),'attempts':1 if l==1 else 100})
rows=[]
for n in range(15,21):
 c=read(f'n{n}-candidate.json');full=read(f'n{n}-search.json');near=read(f'n{n}-near-wronskian.json');lower=read(f'n{n}-lower-multiplicity.json');poles=read(f'n{n}-poles.json');wr=read(f'n{n}-wronskian.json')
 assert near['best']['distance']<float(c['distance'])+1e-8
 assert full['best']['distance']<float(c['distance'])+1e-8
 assert lower['best']['distance']<float(c['distance'])+1e-8
 assert poles['best']['distance']<float(c['distance'])+1e-8
 assert wr['refinedDistinctRoots']==wr['nonzeroDegree']
 # N15 and16 pole reports ARE the original complete logs; N17 was rerun with a high-order-only batch.
 earlier=prior_poles[n]if n==17 else[]
 count=full['attempts']+near['attempts']+lower['attempts']+poles['attempts']+sum(r['attempts']for r in earlier)
 c['searchSummary']={'recordedCompletedStarts':count,'unrestrictedStarts':full['attempts'],'nearWronskianDEStarts':near['attempts'],'lowerMultiplicityDEStarts':lower['attempts'],'poleStarts':poles['attempts']+sum(r['attempts']for r in earlier),'wronskianCandidatesEvaluated':wr['successfulCandidates'],'bestUnrestrictedDistance':full['best']['distance'],'bestNearWronskianDistance':near['best']['distance'],'bestLowerMultiplicityDistance':lower['best']['distance'],'bestPoleDistance':max([poles['best']['distance']]+[r['distance']for r in earlier]),'files':[f'n{n}-search.json',f'n{n}-near-wronskian.json',f'n{n}-lower-multiplicity.json',f'n{n}-poles.json',f'n{n}-wronskian.json'],'additionalCompletedPoleBatch':earlier}
 if n==20:c['searchSummary']['discardedExploratoryBatch']='Eight earlier near-Wronskian DE runs were replaced by eight guarded reruns after a degree-trimming artifact was identified. Those eight preliminary runs are excluded from recordedCompletedStarts; the winning high-precision candidate was unaffected.'
 if n==18:c['searchSummary']['interruptedBatch']='An earlier lower-order pole batch was interrupted before its first100-start group completed; these uncounted partial attempts are excluded.'
 rows.append(c)
 (ROOT/f'n{n}-result.json').write_text(json.dumps(c,indent=2)+'\n')
out={'scope':'One specified maximal1-AC subspace perN; fixed source basis. Maximization candidates, not proven global maxima.','sphereConvention':'(2Re(z),-2Im(z),|z|^2-1)/(1+|z|^2)','results':rows,'recordedCompletedStarts':sum(r['searchSummary']['recordedCompletedStarts']for r in rows),'nonzeroWronskianCandidates':sum(r['searchSummary']['wronskianCandidatesEvaluated']for r in rows),'unrestrictedAndRootConstrainedSearchesFoundNoImprovement':True}
(ROOT/'results.json').write_text(json.dumps(out,indent=2)+'\n')
for r in rows:print(r['n'],r['distance'],r['searchSummary']['recordedCompletedStarts'])
print('total',out['recordedCompletedStarts'])
