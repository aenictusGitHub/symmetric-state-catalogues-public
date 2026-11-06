"""Independently audit the reported AC / equal-weight S^2 design comparison."""
from pathlib import Path
import argparse,csv,hashlib,json
HERE=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--catalogue',type=Path,default=HERE.parents[1]/'data/catalogue.json');parser.add_argument('--baseline',type=Path);args=parser.parse_args()
    data=json.loads(args.catalogue.read_text());benchmark=json.loads((HERE/'benchmark.json').read_text())
    assert benchmark['ambientDimension']==3 and benchmark['sphereDimension']==2
    reference={r['n']:r for r in benchmark['rows']};assert sorted(reference)==list(range(1,101))
    recorded={r['n']:r for r in data['rowMaximumOrders']}
    attached=data['sphericalDesignComparison'];rendered={r['n']:r for r in attached['rows']};assert sorted(rendered)==list(range(1,101))
    rows=[];gaps=[]
    for n in range(1,101):
        entries=[e for e in data['entries'] if e['n']==n]
        exact=max([0]+[e['t'] for e in entries if e['status']=='exists'])
        upper=min(e['t'] for e in entries if e['status']=='impossible')-1
        candidate=max([0]+[s['achievedOrder'] for s in data['states'] if s['n']==n])
        status='proven' if exact==upper else 'conjectured'
        if n==1:assert candidate==exact==upper==0
        else:assert recorded[n]['order']==candidate and recorded[n]['status']==status
        r=reference[n];a=rendered[n];delta=r['order']-candidate
        expected={'n':n,'acOrder':candidate,'acMaximumStatus':status,'acExactAttainmentProven':exact>=candidate,'acRigorousLower':exact,'acRigorousUpper':upper,'designOrder':r['order'],'designConstructionEvidenceCode':r['constructionEvidenceCode'],'designConstructionEvidence':r['constructionEvidence'],'ordersMatch':delta==0,'difference':delta}
        for key,value in expected.items():assert a[key]==value,(n,key,a.get(key),value)
        rows.append(expected)
        if delta:gaps.append(expected)
    assert [r['n'] for r in gaps]==[29,87,94]
    assert attached['mismatchNumbers']==[29,87,94]
    assert all(r['difference']==1 and r['acMaximumStatus']=='conjectured' and r['designConstructionEvidenceCode']=='V3' for r in gaps)
    unchanged=None
    if args.baseline:
        old=json.loads(args.baseline.read_text())
        for key in ['states','entries','proofs','rowMaximumOrders']:
            assert data[key]==old[key],key+' scientific records changed'
        unchanged=True
    report={'passed':True,'checkedDate':'2026-11-06','comparedQubitNumbers':100,'matchingRecordedOrders':97,'mismatchNumbers':[29,87,94],'mismatches':gaps,'comparisonMeaning':'Different published orders; no unequal true maxima are proved by this comparison.','evidenceCaution':'All three AC maxima are conjectured. All three design entries are V3 numerical constructions; no design existence or optimality theorem is inferred.','noDesignBasedACBoundsIntroduced':True,'pureScientificRecordsUnchanged':unchanged,'catalogueSHA256':hashlib.sha256(args.catalogue.read_bytes()).hexdigest(),'benchmarkSHA256':hashlib.sha256((HERE/'benchmark.json').read_bytes()).hexdigest()}
    (HERE/'comparison-audit.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    with (HERE/'comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
