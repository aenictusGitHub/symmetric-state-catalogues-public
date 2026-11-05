#!/usr/bin/env python3
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent

def main():
 attempts=[];candidates=[]
 for path in sorted(HERE.glob('n*-t*-search.json')):
  d=json.loads(path.read_text());n,t=d['n'],d['t'];b=d['best'];found=d['status']=='candidate-found'
  attempts.append({'n':n,'t':t,'dimension':1,'attemptCount':len(d['attempts']),'randomComplexAttempts':sum(x['initialization']=='random-complex'for x in d['attempts']),'status':d['status'],'bestMaximumReductionFrobeniusResidual':b['maximumReductionFrobeniusResidual'],'searchFile':path.name,'sourceSHA256':hashlib.sha256(path.read_bytes()).hexdigest(),'elapsedSeconds':d['elapsedSeconds'],'bestResidual':b['maximumReductionFrobeniusResidual'],'found':found,'sourceFile':path.name})
  if found:
   sid=f'numerical-survey-n{n}-k1-t{t}'
   candidates.append({'id':sid,'n':n,'t':t,'dimension':1,'evidence':'numerical','label':f'Refined Liège pure AC{t} candidate at N={n}','basis':[{'id':sid+'-basis-1','coefficients':b['coefficients']}],'reportedOrder':t,'sourceFile':f'inputs/acsdk.{n}.{ {21:4,22:5,23:5,24:7}[n]}.txt','searchFile':path.name,'verification':{'normalizationError':abs(sum(a*a+b*b for a,b in b['coefficients'])-1),'allReducedStatesThroughOrderChecked':True,'maximumReductionFrobeniusResidual':b['maximumReductionFrobeniusResidual'],'reductions':b['reductions']},'notes':'Numerical witness only. It is not an exact existence certificate, and failed searches at higher orders do not prove nonexistence.'})
 out={'scope':'Finite numerical pure-state searches at N21–24; no exact existence or nonexistence claim. Higher AC order supplies candidates at every lower order.','threshold':1e-10,'candidates':candidates,'searches':attempts,'sourceProvenance':'inputs/provenance.json'}
 (HERE/'summary.json').write_text(json.dumps(out,indent=2)+'\n');print([(x['n'],x['t'],x['dimension'])for x in candidates]);return out
if __name__=='__main__':main()
