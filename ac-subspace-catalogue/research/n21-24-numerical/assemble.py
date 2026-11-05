"""Collect numerical searches; preserve the distinction between evidence and proof."""
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent

def main():
 groups=[];candidates=[];targets=[]
 for label in ['t2','t34','t567']:
  path=HERE/label/'summary.json';d=json.loads(path.read_text());groups.append({'group':label,'source':label+'/summary.json','sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
  for c in d['candidates']:
   assert c['n']in [21,22,23,24]and c['dimension']>=1
   record=dict(c);record['evidence']='numerical';record['sourceGroup']=label;record['groupSummary']=label+'/summary.json'
   record.setdefault('label',f"Numerical search · N={c['n']}, t={c['t']}, dimension {c['dimension']}")
   assert len(record['basis'])==record['dimension'];assert all(len(v['coefficients'])==record['n']+1 for v in record['basis'])
   assert not any('exactCoefficients'in v for v in record['basis'])
   candidates.append(record)
  for target in d.get('targets',d.get('searches',[])):
   r={**target,'sourceGroup':label,'groupSummary':label+'/summary.json'}
   attempts=r.get('attemptCount',r.get('attempts',0));r['attemptCount']=len(attempts)if isinstance(attempts,list)else int(attempts)
   assert isinstance(r['found'],bool)and r['attemptCount']>=0
   source=r.get('sourceFile',r.get('searchFile',label+'/summary.json'))
   r['sourceFile']=source if source.startswith(label+'/') else label+'/'+source
   targets.append(r)
 assert len({c['id']for c in candidates})==len(candidates)
 out={'schemaVersion':1,'date':'2026-11-03','scope':'Numerical whole-subspace searches for N=21–24. Small residuals are numerical evidence, not exact existence certificates; failed searches never prove nonexistence or global optimality.','acceptanceThreshold':1e-10,'groups':groups,'candidates':sorted(candidates,key=lambda c:(c['n'],c['t'],-c['dimension'])),'targets':targets}
 (HERE/'bases.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
 print('Collected',len(candidates),'numerical candidates and',len(targets),'target searches.')
if __name__=='__main__':main()
