#!/usr/bin/env python3
"""Exact multiplicities and checked numerical Majorana rendering of the AC4 line."""
from pathlib import Path
import hashlib,importlib.util,json
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('majorana_helpers',HERE.parents[1]/'scripts/majorana_coordinates.py');helpers=importlib.util.module_from_spec(spec);spec.loader.exec_module(helpers)
def verify():
 raw=(HERE/'basis.json').read_bytes();data=json.loads(raw);base=data['basis'][0];print('building exact geometry',flush=True)
 record=helpers.coordinates(base);assert len(record['points'])==19 and record['partition']==[1]*19
 diag=record['diagnostics'];assert diag['squareFreeFactorizationExact']and diag['zeroMultiplicity']==0 and diag['infinityMultiplicity']==1
 assert diag['projectiveDickeReconstructionError']<1e-12
 output={'schemaVersion':1,'records':{base['id']:record}}
 (HERE/'geometry.json').write_text(json.dumps(output,indent=2)+'\n')
 report={'passed':True,'sourceBasisSha256':hashlib.sha256(raw).hexdigest(),'recordsChecked':1,'allNineteenRootsSimpleExactly':True,'maximumProjectiveDickeReconstructionError':diag['projectiveDickeReconstructionError'],'scope':'Exact square-free polynomial factorization; sphere positions are numerical roots checked by projective Dicke reconstruction.'}
 (HERE/'geometry-verification.json').write_text(json.dumps(report,indent=2)+'\n')
 imported=json.loads((HERE/'pure-import.json').read_text());st=imported['states'][0]
 imported['geometryImports'][st['id']]={'exactCoefficients':st['exactCoefficients'],'points':record['points'],'diagnostics':diag,'sourceSubspaceId':data['id'],'sourceGeometryCertificate':'../ac-subspace-catalogue/research/n19-t4/geometry-verification.json'}
 imported['state']=st;imported['geometry']=imported['geometryImports'][st['id']]
 (HERE/'pure-import.json').write_text(json.dumps(imported,indent=2)+'\n')
 return report
if __name__=='__main__':print(json.dumps(verify(),indent=2))
