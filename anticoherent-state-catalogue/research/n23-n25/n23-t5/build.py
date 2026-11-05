from pathlib import Path
import json,math,argparse
import numpy as np
from scipy.linalg import qr
import mpmath as mp
from fractions import Fraction as Q

here=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description='Discover the rational contraction certificate from the saved numerical N23 source.')
parser.add_argument('--source',type=Path,default=here.parents[2]/'data/catalogue.json',help='AC catalogue JSON; defaults to the current catalogue data.')
args=parser.parse_args()
source=args.source
st=next(s for s in json.loads(source.read_text())['states'] if s['id']=='liege-acs.23.5')
n=23;t=5;d=n+1;m=(t+1)**2
mp.mp.dps=150
scales=[2**((math.comb(n,r).bit_length()-1)//2) for r in range(d)]
z=[mp.mpc(*pair)*scales[r]/mp.sqrt(math.comb(n,r)) for r,pair in enumerate(st['highPrecisionCoefficients'])]
v=[mp.re(x) for x in z]+[mp.im(x) for x in z]
center=[Q(mp.nstr(x,100)) for x in v]
polys=[];labels=[]
for a in range(t+1):
 terms={}
 for l in range(n-t+1):
  r=a+l;w=Q(math.comb(n-t,l),scales[r]**2)
  terms[(r,r)]=w;terms[(r+d,r+d)]=w
 polys.append((terms,-Q(1,(t+1)*math.comb(t,a))));labels.append([a,a,'real'])
for a in range(t+1):
 for b in range(a+1,t+1):
  re={};im={}
  for l in range(n-t+1):
   i=a+l;j=b+l;w=Q(math.comb(n-t,l),scales[i]*scales[j])
   for key,val in [((i,j),w),((i+d,j+d),w)]:re[tuple(sorted(key))]=val
   for key,val in [((i+d,j),w),((i,j+d),-w)]:im[tuple(sorted(key))]=val
  polys.extend([(re,Q(0)),(im,Q(0))]);labels.extend([[a,b,'real'],[a,b,'imag']])
assert len(polys)==m
def ev(p,x):return p[1]+sum(w*x[i]*x[j] for (i,j),w in p[0].items())
def jac(p,x):
 out=[Q(0)]*len(x)
 for (i,j),w in p[0].items():out[i]+=w*x[j];out[j]+=w*x[i]
 return out
F=[ev(p,center) for p in polys]
J=[jac(p,center) for p in polys]
Jn=np.array(J,dtype=float)
_,R,piv=qr(Jn,pivoting=True,mode='economic')
print('rank',np.linalg.matrix_rank(Jn,tol=1e-12),'singular',np.linalg.svd(Jn,compute_uv=False)[-5:],flush=True)
active=[int(i) for i in piv[:m]];fixed=[i for i in range(2*d) if i not in active]
Jm=mp.matrix([[mp.mpf(J[i][j].numerator)/J[i][j].denominator for j in active] for i in range(m)])
inv=Jm**-1
bits=256;den=1<<bits
A=[[Q(int(mp.nint(inv[i,j]*den)),den) for j in range(m)] for i in range(m)]
af=[sum(a*f for a,f in zip(row,F)) for row in A]
eta=max(abs(x) for x in af)
beta=max(sum(abs(Q(i==j)-sum(A[i][k]*J[k][active[j]] for k in range(m))) for j in range(m)) for i in range(m))
anorm=max(sum(abs(a) for a in row) for row in A)
# Sum of the absolute Hessian coefficients, after discarding fixed variables.
aset=set(active)
hnorm=max(sum((2*abs(w) if i==j else 2*abs(w)) for (i,j),w in p[0].items() if i in aset and j in aset) for p in polys)
radius=Q(1,10**65);contraction=beta+anorm*hnorm*radius
print('eta,beta,||A||,||H||,q,r',*[float(x) for x in [eta,beta,anorm,hnorm,contraction,radius]],flush=True)
assert contraction<1 and eta+contraction*radius<radius
cert={'schemaVersion':'rational-quadratic-contraction-v1','n':n,'t':t,'scaleFactors':scales,'variableConvention':'c_r=sqrt(binomial(N,r))/scaleFactors[r]*(x_r+i*x_(N+1+r))','center':[str(x) for x in center],'activeVariables':active,'fixedVariables':fixed,'inverseDyadicBits':bits,'inverseNumerators':[[str(int(x*den)) for x in row] for row in A],'radius':str(radius),'bounds':{k:str(x) for k,x in [('residual',eta),('inverseDefect',beta),('inverseInfinityNorm',anorm),('quadraticHessianBound',hnorm),('contractionFactor',contraction)]},'scope':'Unique exact algebraic zero in the active-variable infinity-norm ball, with twelve other variables fixed to the saved rational center. This certifies existence, not closed-form radical Dicke components.','sourceStateId':st['id']}
(here/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
