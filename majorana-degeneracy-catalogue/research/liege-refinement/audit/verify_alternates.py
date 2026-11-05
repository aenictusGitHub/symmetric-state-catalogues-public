#!/usr/bin/env python3
"""Verify download-only Dicke alternates; deliberately make no partition claim."""
import argparse
import hashlib
import json
from pathlib import Path
import mpmath as mp
from verify_refinement import pairs, direct_reduction, changes, exact_zero, sf


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',type=Path,required=True)
    p.add_argument('--refined',type=Path,required=True)
    p.add_argument('--out',type=Path,default=Path(__file__).with_name('alternate-verification.json'))
    a=p.parse_args()
    mp.mp.dps=110
    original_bytes,refined_bytes=a.original.read_bytes(),a.refined.read_bytes()
    original=json.loads(original_bytes)['records']
    refined=json.loads(refined_bytes)['records']
    by_id={r['id']:r for r in original}
    assert len(by_id)==len(original)
    assert len({r['id'] for r in refined})==len(refined)
    assert set(by_id)=={r['id'] for r in refined}
    results=[]
    for r in refined:
        s=by_id[r['id']]
        assert r['n']==s['n'] and r['reportedOrder']==s['reportedOrder']
        assert r['original90DigitCoefficients']==s['originalCoefficients']
        assert len(r['refinedCoefficients'])==r['n']+1
        assert all(isinstance(x,str) for pair in r['refinedCoefficients'] for x in pair)
        before,after=pairs(s['originalCoefficients']),pairs(r['refinedCoefficients'])
        before_check=direct_reduction(before,r['reportedOrder'])
        after_check=direct_reduction(after,r['reportedOrder'])
        assert mp.mpf(after_check['HSResidual'])<=mp.mpf('1e-60')
        assert mp.mpf(after_check['normError'])<=mp.sqrt(r['reportedOrder']+1)*mp.mpf('1e-60')
        assert r['structuralZeroIndices']==s['structuralZeroIndices']
        assert all(exact_zero(r['refinedCoefficients'][i]) for i in s['structuralZeroIndices'])
        assert r['partition'] is None, 'Alternates must not inherit a coordinate-state partition'
        results.append({'id':r['id'],'sourceObservationId':s['sourceObservationId'],
                        'n':r['n'],'reportedOrder':r['reportedOrder'],'passed':True,
                        'before':before_check,'after':after_check,'change':changes(before,after),
                        'structuralZerosPreservedExactly':True,'partitionStatus':'Not assigned or certified.'})
    summary={'records':len(results),'passedRecords':len(results),
             'maximumBeforeHSResidual':sf(max(mp.mpf(r['before']['HSResidual']) for r in results)),
             'maximumAfterHSResidual':sf(max(mp.mpf(r['after']['HSResidual']) for r in results)),
             'maximumAfterNormError':sf(max(mp.mpf(r['after']['normError']) for r in results)),
             'maximumPhaseAlignedComponentChange':sf(max(mp.mpf(r['change']['maximumPhaseAlignedComponentChange']) for r in results))}
    a.out.write_text(json.dumps({'scientificStatus':'Numerical AC verification only; these unmatched Dicke alternates have no assigned Majorana partition.',
                                'inputs':{'originalFilename':a.original.name,'originalSHA256':hashlib.sha256(original_bytes).hexdigest(),
                                          'refinedFilename':a.refined.name,'refinedSHA256':hashlib.sha256(refined_bytes).hexdigest()},
                                'summary':summary,'records':results},indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
