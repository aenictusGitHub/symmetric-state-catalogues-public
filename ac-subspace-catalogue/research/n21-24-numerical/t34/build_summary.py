#!/usr/bin/env python3
"""Collate completed numerical campaigns without promoting failed searches."""
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent

def main():
 paths=sorted(HERE.glob('n*-t*-k*.json'))+sorted(HERE.glob('structured-n*-t*-k*.json'))+sorted((HERE/'runs').rglob('n*-t*-k*.json'))
 targets={}
 for path in paths:
  data=json.loads(path.read_text());key=(data['n'],data['t'],data['dimension'])
  item=targets.setdefault(key,{'n':key[0],'t':key[1],'dimension':key[2],'k':key[2],'attemptCount':0,'bestResidual':float('inf'),'found':False,'sourceFiles':[],'methods':[]})
  item['attemptCount']+=len(data['attempts']);item['found']|=data['found']
  if data['maxCrossReductionFrobenius']<item['bestResidual']:
   item['bestResidual']=data['maxCrossReductionFrobenius'];item['sourceFile']=str(path.relative_to(HERE))
  item['sourceFiles'].append(str(path.relative_to(HERE)))
  if data['method']not in item['methods']:item['methods'].append(data['method'])
 for item in targets.values():item['interpretation']='A numerical whole-subspace candidate was found; no exact existence claim.'if item['found']else'No candidate at this dimension was found in the recorded finite searches; this is not a nonexistence result.'
 c=json.loads((HERE/'candidates.json').read_text())['candidates'];v=json.loads((HERE/'independent-verification.json').read_text());assert v['passed']
 for rec in c:
  item=targets[rec['n'],rec['t'],rec['dimension']];item['searchBestResidual']=item['bestResidual'];item['bestResidual']=rec['verification']['maxCrossReductionFrobenius'];item['refinedCandidateId']=rec['id'];item['refinementEvaluations']=rec['verification']['nfevRefinement'];item['sourceFile']='candidates.json'
 summary={'scope':'N21–24 AC3 and AC4 numerical dimension searches. AC4 lines are supplied by the complementary pure-state search, not duplicated here.','candidates':c,'targets':list(targets.values()),'attemptCount':sum(x['attemptCount']for x in targets.values()),'independentVerification':'independent-verification.json','highPrecisionVerificationPassed':True,'allSearchesComplete':True,'interpretation':'Numerical dimensions suggest candidates for maxima. Rigorous catalogue intervals and exact maxima remain separate. Every failed larger-dimension search is reported as a failed search, never as dimension zero or an upper-bound proof.'}
 (HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps([{k:x[k]for k in ['n','t','k','attemptCount','bestResidual','found']}for x in targets.values()],indent=2))
if __name__=='__main__':main()
