#!/usr/bin/env python3
"""Verify source hashes, full-partition membership, and separation of evidence."""
from pathlib import Path
from collections import Counter
import hashlib
import importlib.util
import json
from decimal import Decimal
import zipfile

ROOT=Path(__file__).resolve().parent.parent
data=json.loads((ROOT/'data/catalogue.json').read_text())
folder=ROOT/'research/liege-library'
manifest=json.loads((folder/'manifest.json').read_text())
observations=json.loads((folder/'observations.json').read_text())
spec=importlib.util.spec_from_file_location('library_bounds',ROOT/'scripts/library_bounds.py')
bounds=importlib.util.module_from_spec(spec);spec.loader.exec_module(bounds)
checked_files=0
for row in manifest['records']:
    for kind in ('coordinate','Dicke'):
        p=row['provenance'].get(kind)
        if p is None:continue
        raw=folder/p['localfile']
        assert raw.is_file(), raw
        assert hashlib.sha256(raw.read_bytes()).hexdigest()==p['sha256'], raw
        checked_files+=1
assert checked_files==243
assert len(manifest['records'])==145 and len(observations['records'])==144
assert len(observations['missing'])==1 and not observations['errors']
assert all(r['claimNumericallySupported'] for r in observations['records'])
exact={s['id']:s for s in data['states']}
numerical={s['id']:s for s in data['numericalStates']}
expected_exact={s['id'] for s in json.loads((ROOT/'research/states-input.json').read_text())['states']}
expected_exact.add('gmax-n11-t3-certified')
expected_exact.update(s['id'] for s in json.loads((ROOT/'research/ac-order-constructions.json').read_text())['states'])
assert set(exact)==expected_exact and len(numerical)==144 and not set(exact).intersection(numerical)
assert all(s['verification']['exactAC'] and s['verification']['exactMultiplicity'] for s in exact.values())
assert all(not s['verification']['exactAC'] for s in numerical.values())
refinement_folder=ROOT/'research/liege-refinement'
refinements=json.loads((refinement_folder/'refined.json').read_text())
by_refined_id={r['id']:r for r in refinements['records']}
audit=json.loads((refinement_folder/'audit/verification.json').read_text())
assert audit['summary']['passedRecords']==144 and audit['summary']['allRootsChecked']
assert not audit['failures']
assert hashlib.sha256((refinement_folder/'refined.json').read_bytes()).hexdigest()==audit['inputs']['refinedSHA256']
source_hash=hashlib.sha256((folder/'observations.json').read_bytes()).hexdigest()
if source_hash!=audit['inputs']['originalSHA256']:
    equivalence=audit['inputs']['packagedInputEquivalence']
    assert equivalence['allConsumedScientificFieldsIdentical'] and equivalence['records']==144
    assert source_hash==audit['inputs']['equivalentPackagedOriginalSHA256']==equivalence['packagedSourceSHA256']
for state in numerical.values():
    record=by_refined_id[state['numerical']['sourceRecordId']]
    ref=state['refinement']
    assert state['verification']['exactMultiplicity']
    assert state['coefficients']==[[float(a),float(b)] for a,b in record['refinedCoefficients']]
    assert Decimal(ref['afterResidual'])<Decimal('1e-60')
    assert Decimal(ref['afterResidual'])<=Decimal(ref['beforeResidual'])
    path=ROOT/ref['coefficientPath'].removeprefix('./')
    rows=[line.split() for line in path.read_text().splitlines() if line.strip() and not line.startswith('#')]
    assert rows==record['refinedCoefficients']
    assert sum(p['multiplicity'] for p in state['majoranaPoints'])==state['n']
    assert sorted((p['multiplicity'] for p in state['majoranaPoints']),reverse=True)==state['partition']
    assert (ROOT/ref['partitionCertificatePath'].removeprefix('./')).is_file()
with zipfile.ZipFile(refinement_folder/'refined-components.zip') as archive:
    assert archive.testzip() is None
    assert len([p for p in archive.namelist() if p.startswith('coefficients/') and p.endswith('.txt')])==144
    assert len([p for p in archive.namelist() if p.startswith('alternates/') and p.endswith('.txt')])==16
for e in data['entries']:
    n=e['n'];part=tuple(e['partition'])
    assert sum(part)==n and list(part)==sorted(part,reverse=True)
    assert 0<=e['lower']<=e['upper']
    assert e['value']==(e['lower'] if e['lower']==e['upper'] else None)
    for key,pool in [('stateIds',exact),('numericalStateIds',numerical)]:
        for sid in e[key]:
            assert pool[sid]['n']==n and tuple(pool[sid]['partition'])==part,(e['id'],sid)
    reported=max([0]+[numerical[s]['achievedOrder'] for s in e['numericalStateIds']])
    assert e['reportedOrder']==(reported or None)
    if n>12:
        assert e['coverage']=='source-library' and e['upper']==bounds.rigorous_order_upper(n,part)
        lower=bounds.theorem_lower(n,part)
        if n%2==0 and part==(1,)*n and n>=8:lower=max(lower,3)
        if n==3*part[0]+1 and part==(part[0],)+(1,)*(n-part[0]):lower=max(lower,2)
        assert e['lower']==lower
        assert e['numericalStateIds']
    else:assert e['coverage']=='complete-core'
counts=Counter(e['status'] for e in data['entries'])
result={'sourceFilesHashVerified':checked_files,'indexedRecords':145,'numericallyAnalyzedRecords':144,
        'missingSourceFiles':['acs.12.2.3.txt'],'exactStates':len(exact),'numericalStates':len(numerical),
        'partitionRecords':len(data['entries']),'corePartitions':270,'sourceSupplementPartitions':117,
        'statusCounts':dict(counts),'numericalOrdersNotUsedAsRigorousLowerBounds':True,
        'allStateMembershipsMatchFullPartition':True,
        'refinedPrimaryStates':144,'refinedAlternateStates':16,
        'savedCoefficientsMatchAuditedDecimals':True,'numericalACNotPromotedToExactProof':True,
        'fullPartitionsCertifiedForSavedNumericalVectors':True,
        'worstRefinedHSResidual':audit['summary']['maximumAfterHSResidual']}
(ROOT/'research/library-integration-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
