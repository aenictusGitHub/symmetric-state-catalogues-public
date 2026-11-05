#!/usr/bin/env python3
"""Exact attainment of the N18 fourth-order lower bound, reusing the AC5 verifier.
The catalogue coefficients are associated exactly with the checked source ray.
Python standard library only; no numerical tolerance enters the lower-bound proof.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib,importlib.util,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def verify(directory=None):
 here=Path(directory or HERE);source=here.parent/'n18-t5/verify.py'
 spec=importlib.util.spec_from_file_location('n18_exact_ac5_lower',source);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 full=mod.verify();assert full['passed'] and full['exactAchievedOrder']==5
 rho=mod.reduction(4)
 assert len(rho)==5 and all(len(row)==5 for row in rho)
 for a,row in enumerate(rho):
  for b,v in enumerate(row):assert v==({(0,1):F(1,5)} if a==b else {})
 assert mod.distance_squared(rho)==0
 expected=['0']*19
 for r in [0,18]:expected[r]='sqrt(455)/65'
 for r in [5,13]:expected[r]='sqrt(6630)*I/130'
 assert F(455,65**2)==F(14,130) and F(6630,130**2)==F(51,130)
 sub=json.loads((ROOT/'data/catalogue.json').read_text());pure=json.loads((ROOT.parent/'anticoherent-state-catalogue/data/catalogue.json').read_text())
 st=next(x for x in sub['subspaces']if x['id']=='line-acliege-n18-t5');ps=next(x for x in pure['states']if x['id']=='acliege-n18-t5')
 assert (st['n'],st['t'],st['dimension'],st['evidence'])==(18,5,1,'exact')
 assert ps['n']==18 and ps['achievedOrder']==5
 assert st['basis'][0]['exactCoefficients']==ps['exactCoefficients']==expected
 return {'passed':True,'n':18,'t':4,'lowerBound':1,'attainingSubspaceId':st['id'],'pureStateId':ps['id'],'stateLatex':full['stateLatex'],'exactCoefficientSHA256':hashlib.sha256(json.dumps(expected,separators=(',',':')).encode()).hexdigest(),'sourceVerifierSHA256':hashlib.sha256(source.read_bytes()).hexdigest(),'exactNormalization':True,'fourthReducedMatrix':'I_5/5','fourthReducedMatrixEntriesChecked':25,'allReductionsThroughOrderFiveVerified':True,'exactReducedMatrixEntriesChecked':full['exactReducedMatrixEntriesChecked'],'catalogueCoefficientsAssociatedExactly':True,'method':'Reuse of the independent exact normalized-Dicke splitting verifier for the existing AC5 source state, with all 25 fourth-order entries and both catalogue coefficient records checked exactly. Standard library only.','scope':'Attainment of dimension one at fourth order. The upper bound is supplied by the separate auxiliary-spin inertia proof.'}
if __name__=='__main__':
 result=verify();(HERE/'witness-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
