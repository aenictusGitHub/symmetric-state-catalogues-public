#!/usr/bin/env python3
"""Collect successful AC2 numerical witnesses and document failed higher targets."""
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent

def main():
 candidates=[];targets=[]
 audit=json.loads((HERE/'independent-verification.json').read_text());assert audit['passed'];by_n={x['n']:x for x in audit['candidates']}
 for n,k in [(21,4),(22,4),(22,5),(23,4),(23,5)]:
  p=HERE/f'n{n}-k{k}.json';raw=p.read_bytes();d=json.loads(raw);attempts=d['attempts']
  targets.append({'n':n,'t':2,'dimension':k,'attempts':len(attempts),'totalFunctionEvaluations':sum(x['nfev']for x in attempts),'bestResidual':d['bestCrossReductionFrobenius'],'residualDefinition':'Maximum Frobenius norm over every ordered two-qubit cross-reduction deviation, after QR normalization.','found':d['found'],'sourceFile':p.name,'randomSeed':d['randomSeed'],'method':'Full-complex analytic-Jacobian nonlinear least squares; nearby exact-code seeds, independent Gaussian starts, and a long continuation for unsuccessful dimension-five targets.'})
  if d['found']:
   v=by_n[n];assert v['sourceSHA256']==hashlib.sha256(raw).hexdigest()
   candidates.append({'id':d['id'],'n':n,'t':2,'dimension':k,'label':f'Numerical four-dimensional AC2 subspace on {n} qubits','evidence':'numerical','accepted':True,'basis':d['basis'],'references':[{'label':'Reproducible numerical search and coefficients','url':f'research/n21-24-numerical/t2/{p.name}'},{'label':'Independent full-span floating-point verification','url':'research/n21-24-numerical/t2/independent-verification.json'}],'notes':'New numerical full-span AC2 witness, refined to machine precision. Every ordered cross reduction is independently checked. Numerical evidence does not establish exact existence or global optimality.','search':{'sourceFile':p.name,'attempts':len(attempts),'bestSeed':d['bestSeed'],'maxCrossReductionFrobenius':d['bestCrossReductionFrobenius'],'acceptanceThreshold':1e-12},'verification':v})
 out={'orders':[2],'qubitNumbers':[21,22,23,24],'candidates':candidates,'targets':targets,'knownExact':[{'n':24,'t':2,'dimension':5,'id':'algebraic-n24-k5-t2','status':'previously proven maximum','sourceFile':'../../t2-general/n24-k5/basis.json'}],'bestDimensionsFound':{'21':4,'22':4,'23':4,'24':5},'interpretation':'Dimension four is found numerically at N21–23. Dimension five was not found at N22 or N23 in these finite searches; those negative searches are not nonexistence proofs. N24 retains its existing exact optimal dimension-five construction.'}
 (HERE/'summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'candidates':[(x['n'],x['dimension'],x['search']['maxCrossReductionFrobenius'])for x in candidates],'targets':[(x['n'],x['dimension'],x['attempts'],x['totalFunctionEvaluations'],x['found'])for x in targets]},indent=2))
if __name__=='__main__':main()
