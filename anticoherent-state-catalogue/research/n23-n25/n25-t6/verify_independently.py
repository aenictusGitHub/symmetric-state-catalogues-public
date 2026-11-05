#!/usr/bin/env python3
"""Portable exact universal AC6 exclusion check for N=25; standard library only.
Run beside source-identity.json, rational-cuts.json,
verify_source_identity.py, verify_ppt_functional.py and coupled_engine.py.
"""
from pathlib import Path
from importlib.util import spec_from_file_location,module_from_spec
import hashlib,json,sys,time

def load_checker(root,filename,name):
    spec=spec_from_file_location(name,root/filename)
    assert spec and spec.loader
    module=module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return module

def verify(directory=None):
    root=Path(directory or Path(__file__).resolve().parent).resolve()
    sys.path.insert(0,str(root));start=time.monotonic()
    source=load_checker(root,'verify_source_identity.py','n25_source_checker').verify(root)
    assert source['passed'] and source['n']==25 and source['t']==6
    (root/'source-identity-independent.json').write_text(json.dumps(source,indent=2)+'\n')
    ppt=load_checker(root,'verify_ppt_functional.py','n25_ppt_checker').verify(root)
    assert ppt['passed']
    (root/'ppt-independent-verification.json').write_text(json.dumps(ppt,indent=2)+'\n')
    # Both independent reconstructions must concern the exact same payload.
    source_sha=hashlib.sha256((root/'source-identity.json').read_bytes()).hexdigest()
    cuts_sha=hashlib.sha256((root/'rational-cuts.json').read_bytes()).hexdigest()
    assert source['certificateSHA256']==ppt['sourceCertificateSHA256']==source_sha
    assert ppt['rationalCutsSHA256']==cuts_sha
    assert ppt['n']==25 and ppt['copies']==4
    assert ppt['rawSourceFunctionalCoordinatesVerified']==source['variables']==4955
    assert ppt['positiveDyadicRankOneCuts']==81
    files=['source-identity.json','rational-cuts.json']
    report={
        'passed':True,'n':25,'t':6,'copies':4,
        'globalPureStateExclusionProven':True,
        'arbitraryComplexStatesCovered':True,
        'statement':'No pure symmetric 25-qubit state is sixth-order anticoherent. K(25,t)=0 for every t>=6, and the maximum pure AC order is at most five.',
        'attainmentScope':'This verifier checks the upper bound only. Exact fifth-order attainment is certified separately.',
        'method':'Independent exact reconstruction of the physical source trace and all necessary contraction rows; positive source Gram factorizations; independent physical partial-time-reversal functional reconstruction; coefficientwise rational separating identity.',
        'sourceReport':'source-identity-independent.json',
        'pptReport':'ppt-independent-verification.json',
        'pptSummary':{key:ppt[key]for key in ['positiveDyadicRankOneCuts','cutCountsBySplit','parentHighestColumnsReconstructed','parentDescendantNormsChecked','sourceEmbeddingGramEntriesChecked','rawSourceFunctionalCoordinatesVerified']},
        'sourceSummary':{
            'momentDimension':source['inputMomentDimension'],
            'variables':source['variables'],
            'highestWeightVectorsChecked':source['inputHighestWeightVectorsChecked'],
            'physicalCrossGramChecks':source['physicalCrossGramChecks'],
            'contractionRowsReconstructed':source['selectedRowsReconstructed'],
            'positiveSourceBlocks':source['inputPositiveCertificate']['blocks'],
            'positiveRationalLDLPivots':source['inputPositiveCertificate']['positivePivots'],
            'coherentStateControl':source['coherentStateControl']},
        'sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest()for name in files},
        'elapsedSeconds':time.monotonic()-start}
    (root/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    report=verify(sys.argv[1]if len(sys.argv)>1 else None)
    print(json.dumps(report,indent=2))
