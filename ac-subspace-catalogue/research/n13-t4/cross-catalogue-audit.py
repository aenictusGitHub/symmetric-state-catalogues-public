"""Cross-check exact catalogue statuses and numerical local Liege Dicke files.
No floating-point check is promoted to exact existence or nonexistence.
"""
from pathlib import Path
import json,hashlib,re,math
import numpy as np
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[2]
ROOT=PROJECT/'ac-subspace-catalogue'
PURE=PROJECT/'anticoherent-state-catalogue/data/catalogue.json'
SUB=ROOT/'data/catalogue.json'
LIB=Path('/Users/jmartin/Documents/libqiss/ACstates')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def red(c,n,t):
 A=np.array([[c[a+b]*math.sqrt(math.comb(t,a)*math.comb(n-t,b)/math.comb(n,a+b)) for b in range(n-t+1)]for a in range(t+1)])
 return A@A.conj().T

def main():
 pure=json.loads(PURE.read_text());sub=json.loads(SUB.read_text());pm={(e['n'],e['t']):e for e in pure['entries']};sm={(e['n'],e['t']):e for e in sub['entries']}
 contradictions=[];stale=[];stale_sub=[];shared=[]
 for key in sorted(pm.keys()&sm.keys()):
  p,c=pm[key],sm[key];status=p['status'];lower,upper=c['lower'],c['upper']
  if status=='exists' and upper==0 or status=='impossible' and lower>0:contradictions.append({'n':key[0],'t':key[1],'pure':status,'subspaceInterval':[lower,upper]})
  if status in ['unknown','reported-no','reported_no','reportedNo'] and (lower>0 or upper==0):stale.append({'n':key[0],'t':key[1],'pureStatus':status,'newRigorousStatus':'exists' if lower>0 else 'impossible','supportingProofIds':c['proofIds']})
  if status=='impossible' and upper>0:stale_sub.append({'n':key[0],'t':key[1],'subspaceInterval':[lower,upper],'newRigorousUpper':0,'pureProofId':p.get('proofId')})
  if status=='exists' and lower==0:stale_sub.append({'n':key[0],'t':key[1],'subspaceInterval':[lower,upper],'newRigorousLower':1,'pureStateIds':p['stateIds']})
  shared.append({'n':key[0],'t':key[1],'pureStatus':status,'subspaceInterval':[lower,upper]})
 missing=[]
 for n in sub['numbers']:
  if any(key[0]==n for key in pm):continue
  ee=[e for e in sub['entries'] if e['n']==n]
  missing.append({'n':n,'highestExactWitnessOrder':max((e['t'] for e in ee if e['lower']>0),default=0),'provedImpossibleOrders':[e['t'] for e in ee if e['upper']==0],'note':'Size not displayed by the pure-state catalogue.'})
 records=[];numconflicts=[]
 for path in sorted(LIB.glob('acsdk.*.txt'),key=lambda p:tuple(map(int,p.name.split('.')[1:-1]))):
  n,t=map(int,path.name.split('.')[1:3]);a=np.loadtxt(path);assert a.shape==(n+1,2)
  c=a[:,0]+1j*a[:,1];norm=float(np.vdot(c,c).real);c=c/math.sqrt(norm)
  residual=float(np.linalg.norm(red(c,n,t)-np.eye(t+1)/(t+1)))
  nextres=float(np.linalg.norm(red(c,n,t+1)-np.eye(t+2)/(t+2))) if t<n else None
  entry={'filename':path.name,'sha256':sha(path),'n':n,'claimedOrder':t,'normBeforeNormalization':norm,'claimedOrderFrobeniusResidual':residual,'nextOrderFrobeniusResidual':nextres,'numericallySupportedAt1eMinus8':residual<1e-8,'exactExistenceCertificate':False}
  excluded=[q for q in range(1,t+1)if (n,q)in sm and sm[n,q]['upper']==0]
  entry['conflictingProvedExclusions']=excluded
  if excluded and residual<1e-8:numconflicts.append(entry)
  records.append(entry)
 nr=[x for x in records if x['n']==13];assert len(nr)==1 and nr[0]['claimedOrder']==3
 manifest=json.loads((PROJECT/'majorana-degeneracy-catalogue/research/liege-library/manifest.json').read_text())
 deg13=[x for x in manifest['records'] if x['n']==13]
 result={'passed':not contradictions and not numconflicts,'scope':'Consistency audit, not an independent reproof of existing exact certificates.',
 'inputSHA256':{'pureCatalogue':sha(PURE),'subspaceCatalogue':sha(SUB)},'sharedEntriesChecked':len(shared),'contradictions':contradictions,'stalePureStatuses':stale,'staleSubspaceBounds':stale_sub,'sizesMissingFromPureCatalogue':missing,
 'n13SourceStatus':{'DickeFiles':nr,'noFourthOrderDickeFileFound':not (LIB/'acsdk.13.4.txt').exists(),'degenerateFileRecords':deg13,'caution':'In acs.N.t.g.txt, the final index is the reported maximum point multiplicity, not the AC order. Absence of an order-four file is not nonexistence proof.'},
 'localDickeFileCount':len(records),'localNumericalSupportFailures':[x for x in records if not x['numericallySupportedAt1eMinus8']],'localNumericalVsRigorousConflicts':numconflicts,'localNumericalRecords':records,'sharedEntries':shared}
 dest=HERE/'cross-catalogue-audit.json';dest.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['localNumericalRecords','sharedEntries','n13SourceStatus']},indent=2));print('N13',nr)
 return result
if __name__=='__main__':main()
