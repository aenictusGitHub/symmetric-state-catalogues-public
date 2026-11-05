#!/usr/bin/env python3
"""Independent reconstruction check of the refined numerical sphere points."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
from verify_refinement import pairs, changes, sf


def reconstruct(rows):
    q=[mp.mpc(1)]
    norm_errors=[]
    for row in rows:
        x,y,z=(mp.mpf(str(v)) for v in row)
        norm=mp.sqrt(x*x+y*y+z*z)
        norm_errors.append(abs(norm-1))
        x,y,z=x/norm,y/norm,z/norm
        if z>=0:
            a=mp.sqrt((1+z)/2)
            b=(x+mp.j*y)/(2*a)
        else:
            b=mp.sqrt((1-z)/2)
            a=(x-mp.j*y)/(2*b)
        old=q
        q=[mp.mpc(0)]*(len(old)+1)
        for r,coefficient in enumerate(old):
            q[r]+=a*coefficient
            q[r+1]+=b*coefficient
    n=len(rows)
    c=[q[r]/mp.sqrt(math.comb(n,r)) for r in range(n+1)]
    norm=mp.sqrt(mp.fsum(abs(v)**2 for v in c))
    return [v/norm for v in c],max(norm_errors)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--refined',type=Path,required=True)
    p.add_argument('--points',type=Path,required=True)
    p.add_argument('--out',type=Path,default=Path(__file__).with_name('points-verification.json'))
    a=p.parse_args()
    mp.mp.dps=110
    cb,pb=a.refined.read_bytes(),a.points.read_bytes()
    coefficients=json.loads(cb)['records']
    points=json.loads(pb)['records']
    by_id={r['id']:r for r in coefficients}
    assert len({r['id'] for r in points})==len(points)==len(by_id)
    assert set(by_id)=={r['id'] for r in points}
    results=[]
    for r in points:
        source=by_id[r['id']]
        assert r['n']==source['n'] and r['partition']==source['partition']
        assert len(r['highPrecisionRefinedPoints'])==r['n']
        reconstructed,norm_error=reconstruct(r['highPrecisionRefinedPoints'])
        change=changes(pairs(source['refinedCoefficients']),reconstructed)
        assert mp.mpf(change['maximumPhaseAlignedComponentChange'])<mp.mpf('1e-65'), (r['id'],change)
        assert norm_error<mp.mpf('1e-70')
        moved=max(mp.sqrt(mp.fsum((mp.mpf(x)-mp.mpf(y))**2 for x,y in zip(old,new)))
                  for old,new in zip(source['originalPoints'],r['highPrecisionRefinedPoints']))
        results.append({'id':r['id'],'passed':True,'maximumPointNormError':sf(norm_error),
                        'maximumPointDisplacement':sf(moved),
                        'reconstructedCoefficientAgreement':change})
    summary={'records':len(results),'passedRecords':len(results),
             'maximumPointNormError':sf(max(mp.mpf(r['maximumPointNormError']) for r in results)),
             'maximumPointDisplacement':sf(max(mp.mpf(r['maximumPointDisplacement']) for r in results)),
             'maximumReconstructedCoefficientError':sf(max(mp.mpf(r['reconstructedCoefficientAgreement']['maximumPhaseAlignedComponentChange']) for r in results))}
    a.out.write_text(json.dumps({'scientificStatus':'Numerical consistency of displayed/refined sphere points with refined state coefficients; use the separate interval audit for exact multiplicity statements.',
                                'inputs':{'refinedFilename':a.refined.name,'refinedSHA256':hashlib.sha256(cb).hexdigest(),
                                          'pointsFilename':a.points.name,'pointsSHA256':hashlib.sha256(pb).hexdigest()},
                                'summary':summary,'records':results},indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
