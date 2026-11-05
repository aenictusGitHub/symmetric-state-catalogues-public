import json
from pathlib import Path
import mpmath as mp
from analyze_core import reconstruct,rho_residual,clustering
HERE=Path(__file__).parent;mp.mp.dps=90
s=1/mp.sqrt(3);tetra=[(s,s,s),(s,-s,-s),(-s,s,-s),(-s,-s,s)];octa=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
cases=[('tetrahedron',tetra,2,[1]*4),('double tetrahedron',tetra*2,2,[2]*4),('octahedron',octa,3,[1]*6),('double octahedron',octa*2,3,[2]*6),('100-point regular equatorial polygon',[(mp.cos(2*mp.pi*k/100),mp.sin(2*mp.pi*k/100),0) for k in range(100)],1,[1]*100)]
out=[]
for name,points,t,partition in cases:
 points=[tuple(mp.mpf(x) for x in p) for p in points];c,p,e=reconstruct(points);res=rho_residual(c,t);cl=clustering(p);assert cl['partition']==partition and mp.mpf(res['frobeniusError'])<mp.mpf('1e-60');out.append({'test':name,'n':len(points),'t':t,'residual':res,'partitionMatches':True})
(HERE/'analysis-validation.json').write_text(json.dumps({'source':'Locally constructed known constellations; NOT imported source-library data','tests':out},indent=2));print(json.dumps(out,indent=2))
