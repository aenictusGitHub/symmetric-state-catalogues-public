#!/usr/bin/env python3
"""Recompute all full-partition bounds from verified catalogue ingredients.

Portable standard-library check. Place in research/ alongside moment-bounds.json
and partition-audit.json; the default catalogue is ../data/catalogue.json.
This checks combination of theorems/certificates. State reductions, polynomial
multiplicities, and global AC exclusions have their separate exact verifiers.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent

def partitions(total,limit=None):
    if total==0:
        yield ()
        return
    for largest in range(min(total,limit or total),0,-1):
        for tail in partitions(total-largest,largest):
            yield (largest,)+tail

def evaluate(coefficients,x):
    value=F(0)
    for c in reversed(coefficients):value=value*x+c
    return value

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--catalogue',type=Path,default=HERE.parent/'data/catalogue.json')
    parser.add_argument('--moments',type=Path,default=HERE/'moment-bounds.json')
    parser.add_argument('--audit',type=Path,default=HERE/'partition-audit.json')
    parser.add_argument('--output',type=Path,default=HERE/'partition-verification.json')
    args=parser.parse_args()
    catalogue=json.loads(args.catalogue.read_text())
    moment_data=json.loads(args.moments.read_text())
    audit=json.loads(args.audit.read_text())
    moments={(r['n'],r['g'],r['t']):r for r in moment_data['records']}
    for cert in moments.values():
        if cert['status']!='moment_excluded':continue
        n,g,t=cert['n'],cert['g'],cert['t']
        p=list(map(F,cert['polynomialAscending']))
        assert len(p)-1<=t
        assert all(evaluate(p,k)>=0 for k in range(g,n+1))
        average=sum(evaluate(p,k) for k in range(n+1))/F(n+1)
        assert average==F(cert['uniformExpectation'])<0
    global_upper={2:1,3:1,4:2,5:1,6:3,7:2,8:3,9:2,10:3,11:3,12:5}
    legacy=catalogue['legacyEntries']
    for n,upper in global_upper.items():
        rows=[e for e in legacy if e['n']==n]
        assert max(e['t'] for e in rows if e['status']=='proven')==upper
        assert next(e for e in rows if e['t']==upper+1)['status']=='no_ac'
    witnesses=defaultdict(list)
    for state in catalogue['states']:
        lam=tuple(sorted(state['partition'],reverse=True));n=state['n']
        assert sum(lam)==n and all(g>0 for g in lam)
        if n not in global_upper:
            continue  # The larger-size supplement has its own exact-state verifier.
        assert state['achievedOrder']<=global_upper[n]
        witnesses[n,lam].append(state)
    expected={}
    for n in range(2,13):
        for lam in partitions(n):
            g=lam[0];balanced=(n%2==0 and lam==(n//2,n//2))
            first_order=(2*g<n or balanced)
            known=max([0]+[s['achievedOrder'] for s in witnesses[n,lam]])
            if not first_order:
                assert known==0
                lower=upper=0
            else:
                lower=max(1,known);upper=global_upper[n]
                for t in range(1,upper+1):
                    if moments[n,g,t]['status']=='moment_excluded':
                        upper=t-1
                        break
                if n==3*g+1 and lam!=(g,)+(1,)*(n-g):upper=min(upper,1)
            assert lower<=upper
            expected[n,lam]=(lower,upper)
    assert len(expected)==270
    # Larger sizes have a separately marked, source-observed supplement.
    supplied={(e['n'],tuple(e['partition'])):e for e in catalogue['entries'] if e['n']<=12}
    independent={(e['n'],tuple(e['partition'])):e for e in audit['entries']}
    assert len(supplied)==len(independent)==270
    assert supplied.keys()==independent.keys()==expected.keys()
    for key,(lower,upper) in expected.items():
        entry=supplied[key];other=independent[key]
        assert (entry['lower'],entry['upper'])==(other['lower'],other['upper'])==(lower,upper)
        status='no_ac' if upper==0 else 'proven' if lower==upper else 'bounded'
        assert entry['status']==status
        assert entry['value']==(lower if lower==upper else None)
        assert set(entry['stateIds'])=={s['id'] for s in witnesses[key]}
    exact=Counter(lo for lo,hi in expected.values() if lo==hi)
    ranges=Counter(f'{lo}..{hi}' for lo,hi in expected.values() if lo!=hi)
    result={'partitionCount':270,'everyCatalogueBoundMatchesIndependentAudit':True,
            'allExclusionPolynomialsVerifiedRationally':True,
            'exactPartitions':sum(exact.values()),'boundedPartitions':sum(ranges.values()),
            'exactMaximumCounts':dict(sorted(exact.items())),'rangeCounts':dict(ranges),
            'stateCount':sum(st['n']<=12 for st in catalogue['states']),
            'totalExactRepresentativeCount':len(catalogue['states']),
            'scope':'Combines separately verified states and global exclusions with the finite-filter first-order theorem, exact moment polynomials, and second-order saturation partition theorem.'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
