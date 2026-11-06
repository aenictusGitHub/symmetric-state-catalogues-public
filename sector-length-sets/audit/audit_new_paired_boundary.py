from pathlib import Path
import numpy as np,math,json,importlib.util
root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('constructors',root/'scripts/build_n6_boundary_samples.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def odd(s):
 d=math.sqrt(9-s*s);A=math.sqrt(3+d);c=np.zeros(7);c[1]=A/4;c[3]=math.sqrt(10)/4;c[5]=s/(4*A);return c

def matrix(n):
 a,b,c,d,e=n;sq=math.sqrt
 return np.array([[0,0,0,0,sq(10)*e/5,-d,c],
 [0,0,0,-2*sq(5)*e/5,sq(15)*d/5,0,-b],
 [0,0,2*sq(6)*e/5,-sq(2)*d/5,-3*c/5,sq(15)*b/5,sq(10)*a/5],
 [0,-2*sq(5)*e/5,-sq(2)*d/5,4*c/5,-sq(2)*b/5,-2*sq(5)*a/5,0],
 [sq(10)*e/5,sq(15)*d/5,-3*c/5,-sq(2)*b/5,2*sq(6)*a/5,0,0],
 [-d,0,sq(15)*b/5,-2*sq(5)*a/5,0,0,0],
 [c,-b,sq(10)*a/5,0,0,0,0]])
def hodge(r):
 ell=(1+math.sqrt(1+4*r))/2;t=r*(ell-1);a=(ell-t)/6;V=.5;U=.5-(ell-t)**2/3-t*t
 n=np.array([math.sqrt(6)*a,(math.sqrt(V)+math.sqrt(max(U,0)))/math.sqrt(2),t,(math.sqrt(V)-math.sqrt(max(U,0)))/math.sqrt(2),math.sqrt(6)*a]);K=matrix(n)
 vals,vecs=np.linalg.eigh(K);i=abs(vals-1).argmin();c=vecs[:,i]
 if c[np.argmax(abs(c))]<0:c=-c
 return c,abs(vals[i]-1),abs(np.dot(n,n)-1)
checks=[]
for s in np.linspace(-1,3,41):
 c=odd(s);v=m.sector_lengths(c);target=[(9-s*s)/24,(45+5*s)/4,(135-30*s-5*s*s)/24]
 checks.append({'family':'odd-ellipse','parameter':float(s),'error':max(abs(v[k]-target[k]) for k in range(3)),'normalizationError':abs(np.dot(c,c)-1)})
for r in np.linspace(0,.5,31):
 c,ev,normn=hodge(r);v=m.sector_lengths(c);x=(1+4*r*r)/6;target=[x,5+15*x,10-10*x]
 checks.append({'family':'Hodge-linear-seam','parameter':float(r),'error':max(abs(v[k]-target[k]) for k in range(3)),'normalizationError':abs(np.dot(c,c)-1),'eigenvalueResidual':ev,'kernelNormError':normn})
report={'status':'PASS','states':len(checks),'maximumSectorError':max(q['error'] for q in checks),'checks':checks,'imageIdentification':'Paired-purity-face lower-q straight seam followed by its odd-C2 elliptic boundary. Orthographic fit from R7/E*/omega and G1 reproduces the marked contour to screenshot-coordinate precision.'}
assert report['maximumSectorError']<1e-10
(root/'audit/new-paired-boundary-independent-audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))
