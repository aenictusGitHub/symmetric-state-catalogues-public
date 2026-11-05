"""Check attainment, next-order exclusion and subspace consequences at N21/22/25."""
from pathlib import Path
import hashlib,json
here=Path(__file__).resolve().parent;root=here.parents[2]
pure=json.loads((root/'anticoherent-state-catalogue/data/catalogue.json').read_text())
sub=json.loads((root/'ac-subspace-catalogue/data/catalogue.json').read_text());records=[]
for n,t in [(21,4),(22,5),(25,5)]:
    selected=next(e for e in pure['entries']if(e['n'],e['t'])==(n,t));assert selected['status']=='exists'
    state=next(s for s in pure['states']if s['id']==f'radical-n{n}-t{t}')
    assert state['id']in selected['stateIds'] and state['verification']['exactOrderVerified'] and state['verification']['globalMaximumOrderProven']
    maximum=next(m for m in pure['rowMaximumOrders']if m['n']==n)
    assert maximum['order']==t and maximum['exactAttainmentProven'] and maximum['upperBoundProven'] and maximum['status']=='proven'
    for e in pure['entries']:
        if e['n']==n and e['t']>t:assert e['status']=='impossible'
    for e in sub['entries']:
        if e['n']!=n:continue
        if e['t']<=t:assert e['lower']>=1
        if e['t']>t:assert e['lower']==e['upper']==0
    for s in pure['states']:
        if s['n']==n and s.get('status')=='numerical':assert not s['verification']['globalMaximumOrderProven']
    records.append({'n':n,'maximumACOrder':t,'provedMaximum':True,'exactAttainmentAndNextOrderExclusion':True})
checks=[('n21-n22/n21-t5','certificate.json','certificateMetadataSHA256'),('n21-n22/n22-t6','certificate.json','certificateSHA256')]
for folder,filename,key in checks:
    d=here.parent/folder;report=json.loads((d/'independent-verification.json').read_text());assert report['passed'];assert report[key]==hashlib.sha256((d/filename).read_bytes()).hexdigest()
d=here.parent/'n23-n25/n25-t6';report=json.loads((d/'independent-verification.json').read_text());assert report['passed']and report['globalPureStateExclusionProven']and report['arbitraryComplexStatesCovered']
for fn,sha in report['sha256'].items():assert hashlib.sha256((d/fn).read_bytes()).hexdigest()==sha
line=next(e for e in sub['entries']if(e['n'],e['t'])==(22,5));assert line['lower']==line['upper']==1
report={'passed':True,'records':records,'subspaceConsequences':['K(21,t)=0 for t>=5','K(25,t)=0 for t>=6'],'additionalExactMaximum':{'n':22,'t':5,'dimension':1}}
(here/'catalogue-consistency.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
