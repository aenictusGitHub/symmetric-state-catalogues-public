"""Independent reconstruction of every represented state from spherical stars.
The audit reads catalogue coefficients, never the coordinate generator internals.
"""
from pathlib import Path
import json, math, hashlib
import mpmath as mp
mp.mp.dps=80
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT.parent

def reconstruct(n,points):
    p=[mp.mpc(1)]
    for star in points:
        x,y,z=(mp.mpf(str(star[k])) for k in ('x','y','z'))
        v=mp.sqrt(x*x+y*y+z*z);x/=v;y/=v;z/=v
        if z<=0:
            b=mp.sqrt((1-z)/2);a=(x-1j*y)/(2*b)
        else:
            a=mp.sqrt((1+z)/2);b=(x+1j*y)/(2*a)
        for j in range(star.get('multiplicity',1)):
            nxt=[mp.mpc(0) for _ in range(len(p)+1)]
            for r,c in enumerate(p):nxt[r]-=a*c;nxt[r+1]+=b*c
            p=nxt
    assert len(p)==n+1
    c=[(-1)**r*p[r]/mp.sqrt(mp.binomial(n,r)) for r in range(n+1)]
    norm=mp.sqrt(sum(abs(z)**2 for z in c));return [z/norm for z in c]

def normalized(c):
    c=[mp.mpc(*z) if isinstance(z,list) else mp.mpc(z) for z in c]
    norm=mp.sqrt(sum(abs(z)**2 for z in c));return [z/norm for z in c]

def residual(a,b):
    overlap=sum(mp.conj(x)*y for x,y in zip(a,b));phase=overlap/abs(overlap) if overlap else 1
    return float(mp.sqrt(sum(abs(y-phase*x)**2 for x,y in zip(a,b))))

reports=[]
for name,key in [(ROOT.name, 'states' if 'states' in json.loads((ROOT/'data/catalogue.json').read_text()) else 'entries')]:
    root=BASE/name
    science=json.loads((root/'data/catalogue.json').read_text())
    derived=json.loads((root/'data/majorana.json').read_text())
    records=derived['records']; records=records if isinstance(records,dict) else {r['id']:r for r in records}
    checks=[];missing=[]
    for entry in science[key]:
        id=entry['id'];n=entry['n'];rec=records.get(id)
        if not rec:missing.append(id);continue
        points=rec['points']
        assert rec['n']==n,(id,'wrong n')
        assert sum(p['multiplicity'] for p in points)==n,(id,'wrong weight')
        assert all(isinstance(p['multiplicity'],int) and p['multiplicity']>0 for p in points),(id,'invalid multiplicity')
        norms=[abs(math.sqrt(sum(p[k]**2 for k in ('x','y','z')))-1) for p in points]
        assert max(norms)<1e-12,(id,'unit sphere')
        c=reconstruct(n,points);target=normalized(entry['coefficients']);err=residual(target,c)
        assert err<1e-9,(id,err)
        checks.append({'id':id,'n':n,'distinctPositions':len(points),'partition':sorted([p['multiplicity'] for p in points],reverse=True),'phaseAlignedVectorError':err,'unitSphereError':max(norms)})
    assert not missing,(name,missing)
    reports.append({'catalogue':name,'total':len(science[key]),'checked':len(checks),'missing':missing,'maxPhaseAlignedVectorError':max(r['phaseAlignedVectorError'] for r in checks),'maxUnitSphereError':max(r['unitSphereError'] for r in checks),'catalogueSha256':hashlib.sha256((root/'data/catalogue.json').read_bytes()).hexdigest(),'majoranaSha256':hashlib.sha256((root/'data/majorana.json').read_bytes()).hexdigest(),'records':checks})
report={'method':'80-digit homogeneous polynomial reconstruction from unit-sphere points, with exact projective infinity chart, repeated factors weighted by multiplicity, binomial Dicke conversion, normalization and phase alignment against the original stored Dicke vector. Independent of coordinate generator.','passed':True,'threshold':1e-9,'catalogues':reports}
Path(__file__).with_name('coordinates-audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({**report,'catalogues':[{k:v for k,v in r.items() if k!='records'} for r in reports]},indent=2))
