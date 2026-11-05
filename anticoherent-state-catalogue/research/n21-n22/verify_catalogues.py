"""Check that both catalogue snapshots use the verified N21/N22 results."""
from pathlib import Path
import hashlib,json

here=Path(__file__).resolve().parent
root=here.parents[2]
pure=json.loads((root/'anticoherent-state-catalogue/data/catalogue.json').read_text())
sub=json.loads((root/'ac-subspace-catalogue/data/catalogue.json').read_text())
records=[]
for n,lower,upper in [(21,4,5),(22,5,5)]:
    selected=next(e for e in pure['entries'] if (e['n'],e['t'])==(n,lower))
    assert selected['status']=='exists'
    state=next(s for s in pure['states'] if s['id']==f'radical-n{n}-t{lower}')
    assert state['id'] in selected['stateIds']
    assert state['verification']['exactOrderVerified']
    assert state['verification']['globalMaximumOrderProven']==(n==22)
    maximum=next(m for m in pure['rowMaximumOrders'] if m['n']==n)
    assert maximum['exactAttainmentProven']
    assert maximum['status']==('proven' if n==22 else 'conjectured')
    for e in pure['entries']:
        if e['n']==n and e['t']>upper:assert e['status']=='impossible'
    for e in sub['entries']:
        if e['n']!=n:continue
        if e['t']<=lower:assert e['lower']>=1
        if e['t']>upper:assert e['lower']==e['upper']==0
    for directory in ([f'n{n}-t6']):
        cert=here/directory/'certificate.json'
        report=json.loads((cert.parent/'independent-verification.json').read_text())
        assert report['passed'] and report['arbitraryComplexStatesCovered']
        assert report['certificateSHA256']==hashlib.sha256(cert.read_bytes()).hexdigest()
    records.append({'n':n,'maximumACOrderLower':lower,'maximumACOrderUpper':upper,'provedMaximum':n==22})
line=next(e for e in sub['entries'] if (e['n'],e['t'])==(22,5))
assert line['lower']==line['upper']==1
gap=next(e for e in pure['entries'] if (e['n'],e['t'])==(21,5))
assert gap['status']=='unknown'
report={'passed':True,'records':records,'subspaceConsequence':{'n':22,'t':5,'exactMaximumDimension':1},'remainingPureStateQuestion':{'n':21,'t':5,'status':'unresolved'}}
(here/'catalogue-consistency.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
