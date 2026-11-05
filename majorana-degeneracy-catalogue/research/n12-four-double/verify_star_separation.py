#!/usr/bin/env python3
"""Interval-certify distinct stars for serialized geometric-search parameters.

The exact state meant here is the normalized symmetrized product of the eight
spinors (cos(theta/2),exp(i phi)sin(theta/2)), at the stored decimal angles,
with multiplicities(2,2,2,2,1,1,1,1). This certifies its full partition only;
the separate AC residual is numerical and does not prove exact anticoherence.
"""
import argparse
import hashlib
import itertools
import json
from decimal import Decimal,localcontext,ROUND_FLOOR
from pathlib import Path
import mpmath as mp

HERE=Path(__file__).resolve().parent


def floor_string(value):
    sign,mantissa,exponent,_=value._mpf_
    numerator=(-1 if sign else 1)*mantissa*(1<<max(exponent,0))
    denominator=1<<max(-exponent,0)
    with localcontext() as context:
        context.prec=35;context.rounding=ROUND_FLOOR
        return str(Decimal(numerator)/Decimal(denominator))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputs',nargs='*',type=Path)
    p.add_argument('--output',type=Path,default=HERE/'verified-star-partitions.json');a=p.parse_args()
    files=a.inputs or [HERE/'geometric-search-results.json',HERE/'geometric-double-bound-results.json']
    mp.mp.dps=90;mp.iv.dps=90;records=[];sources=[]
    mapping=[(1,0)]+[(i,j) for i in range(2,8) for j in (0,1)]
    for path in files:
        raw=path.read_bytes();data=json.loads(raw)
        sources.append({'filename':path.name,'sha256':hashlib.sha256(raw).hexdigest(),'runs':len(data['runs'])})
        for index,run in enumerate(data['runs']):
            angles=[[mp.iv.mpf(0),mp.iv.mpf(0)] for _ in range(8)]
            for value,(i,j) in zip(run['variables'],mapping):angles[i][j]=mp.iv.mpf(str(value))
            stars=[[mp.iv.sin(th)*mp.iv.cos(ph),mp.iv.sin(th)*mp.iv.sin(ph),mp.iv.cos(th)] for th,ph in angles]
            distances=[]
            for i,j in itertools.combinations(range(8),2):
                squared=sum((a-b)**2 for a,b in zip(stars[i],stars[j]))
                lower=mp.mpf(mp.iv.sqrt(squared)._mpi_[0])
                assert lower>0,(path,index,i,j)
                distances.append(lower)
            records.append({'source':path.name,'runIndex':index,'seed':run['seed'],
                            'minimumChordDistanceCertifiedLower':floor_string(min(distances)),
                            'fullPartition':[2,2,2,2,1,1,1,1],'partitionCertified':True,
                            'ACStatus':'Numerical residual only; exact3AC not established.'})
    result={'scientificStatus':'Exact full multiplicity partition for the states defined by the serialized angles, checked with outward interval trigonometry. No exact AC claim.',
            'precisionDigits':90,'sources':sources,'summary':{'states':len(records),'allPartitionsCertified':True,
            'minimumCertifiedChordDistance':str(min(Decimal(r['minimumChordDistanceCertifiedLower']) for r in records))},'records':records}
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(result['summary'])


if __name__=='__main__':main()
