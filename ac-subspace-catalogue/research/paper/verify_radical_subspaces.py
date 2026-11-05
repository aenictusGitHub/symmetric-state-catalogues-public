#!/usr/bin/env python3
"""Independent full-subspace verification of the radical reconstructions.
Usage: python3 verify_radical_subspaces.py [input.json]
The verification imports only the generic Dicke splitting check, not the derivation.
"""
from pathlib import Path
import json,hashlib,sys
from verify_subspaces import symbolic_check
ROOT=Path(__file__).resolve().parent
source=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'radical-subspaces-input.json'
raw=source.read_bytes();data=json.loads(raw);reports=[]
for rec in data['subspaces']:
 rows=[v['exactCoefficients'] for v in rec['basis']]
 result=symbolic_check(rows,rec['n'],rec['order'])
 assert result['passed'] and result['orthonormal'] and result['fullRank'],(rec['id'],result)
 achieved=rec['order'];higher=symbolic_check(rows,rec['n'],achieved+1)
 while higher['passed'] and achieved<rec['n']//2:
  achieved+=1;higher=symbolic_check(rows,rec['n'],achieved+1)
 bound_paper=(rec['n']-achieved+1)//(achieved+1)
 bound_tensor=(rec['n']+1)//((achieved//2+1)**2)
 report={'id':rec['id'],'n':rec['n'],'dimension':rec['dimension'],'claimedOrder':rec['order'],'verifiedAchievedOrder':achieved,**result,'firstFailedOrder':achieved+1,'failureAtNextOrder':higher['violations'][0], 'rigorousMaximumDimensionUpperBound':min(bound_paper,bound_tensor),'maximalDimensionProven':rec['dimension']==min(bound_paper,bound_tensor)}
 if report['maximalDimensionProven']:
  report['maximalityProof']=f"The exact basis spans a {rec['dimension']}-dimensional {achieved}-AC subspace. Eq. (52) bounds any such subspace by floor(({rec['n']}−{achieved}+1)/({achieved}+1))={bound_paper}; the tensor-Gram bound is {bound_tensor}. The attained dimension equals the best bound."
 reports.append(report)
output={'passed':True,'sourceSha256':hashlib.sha256(raw).hexdigest(),'method':'Independent normalized-Dicke partial-trace splitting; all ordered basis cross matrices and all their scalar entries are checked exactly, as are basis orthogonality and rank.','subspaces':reports}
(ROOT/'radical-independent-verification.json').write_text(json.dumps(output,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(output,indent=2,ensure_ascii=False))
