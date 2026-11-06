"""Import the published Hardin--Sloane exact-N spherical-design benchmark.

Only factual N/order/evidence fields are retained. Construction evidence is
separate from optimality: the source describes its strengths as conjectured.
No spherical-design record becomes an anticoherence upper bound.
"""
from pathlib import Path
from html.parser import HTMLParser
import argparse,hashlib,json,re,csv
HERE=Path(__file__).resolve().parent
URL='https://neilsloane.com/sphdesigns/index.html'
class Preformatted(HTMLParser):
    def __init__(self):super().__init__();self.depth=0;self.blocks=[];self.current=[]
    def handle_starttag(self,tag,attrs):
        if tag=='pre':
            if self.depth==0:self.current=[]
            self.depth+=1
    def handle_endtag(self,tag):
        if tag=='pre' and self.depth:
            self.depth-=1
            if not self.depth:self.blocks.append(''.join(self.current))
    def handle_data(self,data):
        if self.depth:self.current.append(data)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--html',type=Path,required=True);args=parser.parse_args()
    raw=args.html.read_bytes();p=Preformatted();p.feed(raw.decode('utf-8'))
    table=next(block for block in p.blocks if 'tau(N)' in block)
    rows=[]
    evidence={'V1':'algebraic-existence','V2':'validated-existence','V3':'numerical-construction'}
    for line in table.splitlines():
        match=re.match(r'^\s*(\d+)\s+(\d+)\s+(.+?)\s*$',line)
        if not match:continue
        n,order=int(match[1]),int(match[2]);tail=match[3]
        if not 1<=n<=100:continue
        code=tail.split()[0] if tail.split()[0] in evidence else 'reference'
        reference=(re.match(r'\[([^]]+)\]',tail).group(1) if tail.startswith('[') else 'McLaren 1963' if tail.startswith('McL63') else None)
        row={'n':n,'order':order,'constructionEvidenceCode':code,'constructionEvidence':evidence.get(code,'literature-construction-cited'),'optimality':'putative; not inferred from construction evidence','coordinateUrl':f'https://neilsloane.com/sphdesigns/dim3/des.3.{n}.{order}'}
        if reference:row['constructionReference']=reference
        rows.append(row)
    assert [r['n'] for r in rows]==list(range(1,101))
    assert [(r['n'],r['order'],r['constructionEvidenceCode']) for r in rows if r['n'] in [29,87,94]]==[(29,6,'V3'),(87,12,'V3'),(94,13,'V3')]
    data={'schemaVersion':1,'retrievedDate':'2026-11-06','ambientDimension':3,'sphereDimension':2,'pointConvention':'Exactly N distinct points on S^2 with equal weights.','source':{'title':'Hardin–Sloane: N-point spherical designs in three dimensions','authors':['R. H. Hardin','N. J. A. Sloane'],'url':URL,'paperUrl':'https://doi.org/10.1007/BF02711518','tableLastModified':'2002-10-02','retrievedDate':'2026-11-06','sourceHTMLSHA256':hashlib.sha256(raw).hexdigest(),'status':'Published putatively optimal strengths, not a complete table of proved maxima.','evidenceCodes':{'V1':'Algebraic existence proof of the construction.','V2':'Interval existence proof of the construction.','V3':'Numerical construction with source-reported discrepancy at most 1e-26.','reference':'Existence construction credited to the cited literature; not independently re-proved here.'}},'scope':'A reference comparison of reported orders. These strengths are not used as AC existence proofs or AC upper bounds. Construction evidence does not imply optimality.','rows':rows}
    (HERE/'benchmark.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    with (HERE/'benchmark.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['n','order','constructionEvidenceCode','constructionEvidence']);w.writeheader();w.writerows({k:r[k] for k in w.fieldnames} for r in rows)
    print('Imported100 exact-N benchmark rows; all3 reported-order gaps useV3 design constructions.')
if __name__=='__main__':main()
