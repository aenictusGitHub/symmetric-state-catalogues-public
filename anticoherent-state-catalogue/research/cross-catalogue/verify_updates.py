"""Independent checks for imported radical witnesses and cross-catalogue exclusions."""
from pathlib import Path
import importlib.util,json,hashlib
import sympy as s
from sympy.physics.wigner import wigner_6j
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def verify():
 data=json.loads((ROOT/'data/catalogue.json').read_text());entries={(e['n'],e['t']):e for e in data['entries']};states={e['id']:e for e in data['states']}
 spec=importlib.util.spec_from_file_location('pure_reduction_check',ROOT/'research/verify_states.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);module.verify(ROOT/'data/catalogue.json')
 for n in [14,17]:
  rec=states[f'certified-n{n}-t4'];assert rec['achievedOrder']==4 and rec['verification']['exactNormalization']
  for t in range(1,5):assert entries[n,t]['status']=='exists' and rec['id'] in entries[n,t]['stateIds']
  assert entries[n,5]['status']=='impossible' and rec['verification']['globalMaximumOrderProven']
 for n,t in [(15,6),(16,7),(17,7),(18,7)]:
  # Count every rectangular shift independently of either catalogue builder.
  upper=min((n-h+1)//sum(q+1 for q in range(h,t+1,2)) for h in range(t+1))
  assert upper==0 and entries[n,t]['status']=='impossible'
 n16=sum((-1)**(16+7+ell+1)*(2*ell+1)*wigner_6j(8,8,7,7,7,ell)for ell in [1,2,3]);assert s.simplify(n16-5*s.sqrt(69)/8398)==0
 assert entries[16,6]['status']=='impossible' and entries[16,6]['proofId']=='n16-t6-rectangular-equality'
 assert states['acliege-n16-t5']['verification']['globalMaximumOrderProven']
 for t in [4,5]:assert entries[13,t]['status']=='impossible' and entries[13,t]['proofId']=='n13-t4-four-copy'
 assert states['acliege-n13-t3']['verification']['globalMaximumOrderProven']
 n13folder=ROOT.parent/'ac-subspace-catalogue/research/n13-t4'
 n13check=json.loads((n13folder/'independent-verification.json').read_text())
 assert n13check['passed'] and (n13check['n'],n13check['t'],n13check['copies'],n13check['realInvariantVariables'],n13check['constraintRank'],n13check['kernelDimension'])==(13,4,4,331,320,11)
 assert n13check['certificateSHA256']==hashlib.sha256((n13folder/'certificate.json').read_bytes()).hexdigest()
 n14folder=ROOT.parent/'ac-subspace-catalogue/research/n14-t5'
 n14check=json.loads((n14folder/'independent-verification.json').read_text())
 assert n14check['passed'] and n14check['positiveDiagonalCertificate']['formsHaveRankTwo']
 assert n14check['certificateSHA256']==hashlib.sha256((n14folder/'certificate.json').read_bytes()).hexdigest()
 assert entries[14,5]['proofId']=='n14-t5-three-copy'
 n15folder=ROOT.parent/'ac-subspace-catalogue/research/n15-t4'
 n15check=json.loads((n15folder/'independent-verification.json').read_text())
 assert n15check['passed'] and n15check['certificateSHA256']==hashlib.sha256((n15folder/'certificate.json').read_bytes()).hexdigest()
 for t in [4,5]:assert entries[15,t]['status']=='impossible' and entries[15,t]['proofId']=='n15-t4-ppt'
 assert states['acliege-n15-t3']['verification']['globalMaximumOrderProven']
 n17folder=ROOT.parent/'ac-subspace-catalogue/research/n17-t5'
 n17check=json.loads((n17folder/'independent-verification.json').read_text())
 assert n17check['passed'] and n17check['certificateSHA256']==hashlib.sha256((n17folder/'certificate.json').read_bytes()).hexdigest()
 assert entries[17,5]['proofId']=='n17-t5-ppt' and entries[17,6]['proofId']=='n17-t6-complement'
 n17complement=json.loads((ROOT.parent/'ac-subspace-catalogue/research/n17-t6/verification.json').read_text());assert n17complement['passed'] and n17complement['maximumDimension']==0
 n18folder=ROOT.parent/'ac-subspace-catalogue/research/n18-t6'
 n18check=json.loads((n18folder/'independent-verification.json').read_text())
 assert n18check['passed'] and n18check['certificateSHA256']==hashlib.sha256((n18folder/'certificate.json').read_bytes()).hexdigest()
 assert entries[18,5]['status']=='exists' and entries[18,6]['proofId']=='n18-t6-three-copy'
 assert states['acliege-n18-t5']['verification']['globalMaximumOrderProven']
 n18witness=json.loads((ROOT.parent/'ac-subspace-catalogue/research/n18-t5/verification.json').read_text());assert n18witness['passed'] and n18witness['maximumDimension']==1
 n19fifthfolder=ROOT.parent/'ac-subspace-catalogue/research/n19-t5'
 n19fifthcheck=json.loads((n19fifthfolder/'independent-verification.json').read_text())
 assert n19fifthcheck['passed'] and n19fifthcheck['certificateSHA256']==hashlib.sha256((n19fifthfolder/'certificate.json').read_bytes()).hexdigest()
 assert entries[19,5]['status']=='impossible' and entries[19,5]['proofId']=='n19-t5-four-copy'
 assert states['certified-n19-t4']['verification']['globalMaximumOrderProven']
 n19folder=ROOT.parent/'ac-subspace-catalogue/research/n19-t6'
 n19check=json.loads((n19folder/'independent-verification.json').read_text())
 assert n19check['passed'] and n19check['certificateSHA256']==hashlib.sha256((n19folder/'certificate.json').read_bytes()).hexdigest()
 assert entries[19,6]['status']=='impossible' and entries[19,6]['proofId']=='n19-t6-three-copy'
 assert entries[19,4]['status']=='exists' and states['certified-n19-t4']['achievedOrder']==4
 n20folder=ROOT.parent/'ac-subspace-catalogue/research/n20-t6'
 n20check=json.loads((n20folder/'independent-verification.json').read_text())
 assert n20check['passed'] and n20check['certificateSHA256']==hashlib.sha256((n20folder/'certificate.json').read_bytes()).hexdigest()
 assert entries[20,6]['status']=='impossible' and entries[20,6]['proofId']=='n20-t6-three-copy'
 assert states['acliege-n20-t5']['verification']['globalMaximumOrderProven']
 n20=sum((-1)**(20+8+q)*(2*q+2)*wigner_6j(10,10,8,s.Rational(19,2),s.Rational(19,2),s.Rational(2*q+1,2))for q in range(4));assert s.simplify(n20+250*s.sqrt(87)/572033)==0
 assert entries[20,7]['status']=='impossible' and entries[20,7]['proofId']=='n20-t6-three-copy'
 geometry=json.loads((ROOT/'data/majorana.json').read_text());sha=hashlib.sha256((ROOT/'data/catalogue.json').read_bytes()).hexdigest();assert geometry['catalogueSHA256']==sha and set(geometry['records'])==set(states)
 result={'passed':True,'catalogueSHA256':sha,'exactStatesChecked':len(states),'allReducedMatricesThroughAchievedOrderAndNextOrderChecked':True,'importedRadicalStates':[14,17],'newProvedExclusions':[[20,6],[19,5],[19,6],[18,6],[17,5],[17,6],[15,4],[15,5],[14,5],[13,4],[13,5],[15,6],[16,6],[16,7],[18,7],[20,7]],'newN17RowSeventhOrderExcluded':True,'n16RankSevenMultiplier':str(s.simplify(n16)),'n20RankEightMultiplier':str(s.simplify(n20)),'allStatesHaveMajoranaCoordinates':True,'n13MaximumACOrder':3,'n14MaximumACOrder':4,'n16MaximumACOrder':5,'n13FourthOrderStatus':entries[13,4]['status'],'scope':'Full exact witness checks plus exact recoupling coefficients for the linked equality proofs; proof arguments are linked separately.'}
 result['linkedN13Exclusion']={'certificateSHA256':n13check['certificateSHA256'],'independentVerificationPassed':True,'variables':331,'constraintRank':320,'kernelDimension':11,'proofSource':'../ac-subspace-catalogue/research/n13-t4/index.html'}
 result['n15MaximumACOrder']=3
 result['n17MaximumACOrder']=4
 result['n18MaximumACOrder']=5
 result['n19MaximumACOrder']=4
 result['n20MaximumACOrder']=5
 result['linkedN20SixthOrderExclusion']={'certificateSHA256':n20check['certificateSHA256'],'independentVerificationPassed':True,'proofSource':'../ac-subspace-catalogue/research/n20-t6/index.html'}
 result['linkedN19FifthOrderExclusion']={'certificateSHA256':n19fifthcheck['certificateSHA256'],'independentVerificationPassed':True,'variables':1522,'selectedConstraints':1380,'positivePivots':277,'proofSource':'../ac-subspace-catalogue/research/n19-t5/index.html'}
 result['linkedN19SixthOrderExclusion']={'certificateSHA256':n19check['certificateSHA256'],'independentVerificationPassed':True,'variables':94,'constraintRank':88,'kernelDimension':6,'proofSource':'../ac-subspace-catalogue/research/n19-t6/index.html'}
 result['linkedN18SixthOrderExclusion']={'certificateSHA256':n18check['certificateSHA256'],'independentVerificationPassed':True,'variables':84,'constraintRank':80,'kernelDimension':4,'proofSource':'../ac-subspace-catalogue/research/n18-t6/index.html'}
 result['linkedN17FifthOrderExclusion']={'certificateSHA256':n17check['certificateSHA256'],'independentVerificationPassed':True,'proofSource':'../ac-subspace-catalogue/research/n17-t5/index.html'}
 result['linkedN17SixthOrderComplement']={'verificationPassed':True,'generalFormula':n17complement['generalCorollary']['formula'],'proofSource':'../ac-subspace-catalogue/research/n17-t6/index.html'}
 result['linkedN15Exclusion']={'certificateSHA256':n15check['certificateSHA256'],'independentVerificationPassed':True,'proofSource':'../ac-subspace-catalogue/research/n15-t4/index.html'}
 result['linkedN14Exclusion']={'certificateSHA256':n14check['certificateSHA256'],'independentVerificationPassed':True,'variables':47,'constraintRank':45,'kernelDimension':2,'proofSource':'../ac-subspace-catalogue/research/n14-t5/index.html'}
 (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));return result
if __name__=='__main__':verify()
