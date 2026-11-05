"""Independent full degree-eight SOS verification over Q(sqrt(61)).
Only Python's standard library; exact polynomial expansion and matrix tests.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from math import comb,lcm,isqrt
from functools import lru_cache
import json,hashlib,sys
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
BASE=Path(sys.argv[1])if len(sys.argv)>1 else Path(__file__).resolve().parent
basis_file=BASE/'basis.json';cert_file=BASE/'sos-certificate.json'
basis=json.loads(basis_file.read_text());cert=json.loads(cert_file.read_text())
Z=(F(0),F(0));ONE=(F(1),F(0))
def read(v):return tuple(map(F,v))
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def neg(x):return (-x[0],-x[1])
def mul(x,y):return (x[0]*y[0]+61*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def div(x,y):
 d=y[0]**2-61*y[1]**2;assert d
 return mul(x,(y[0]/d,-y[1]/d))
def scale(x,q):return (x[0]*q,x[1]*q)
def positive(x):
 a,b=x
 if not b:return a>0
 if not a:return b>0
 if a*b>0:return a>0
 return a*a>61*b*b if a>0 else a*a<61*b*b
assert mul((F(0),F(1)),(F(0),F(1)))==(F(61),F(0))
bound=(F(130928,1461),F(6832,1461));assert read(cert['bound'])==bound and cert['n']==9 and cert['sector']==5
# Base-nine exponent encoding: every monomial has degree at most eight,
# so integer addition exactly multiplies monomials without digit carries.
def encode(mon):return sum(9**i for i in mon)
SPLIT=9**10
def conjugate(m):return m//SPLIT+(m%SPLIT)*SPLIT
def product(p,q):
 out=defaultdict(int)
 for m,c in p.items():
  for n,d in q.items():out[m+n]+=c*d
 return {m:c for m,c in out.items()if c}
def purity(t):
 p=defaultdict(int)
 for i in range(t+1):
  for j in range(t+1):
   for a in range(10-t):
    for b in range(10-t):p[encode((i+a,10+j+a,j+b,10+i+b))]+=comb(t,i)*comb(t,j)*comb(9-t,a)*comb(9-t,b)
 return p
Q={encode((i,10+i)):comb(9,i)for i in range(10)};Q2=product(Q,Q);Q4=product(Q2,Q2);P5=defaultdict(int)
for t in range(6):
 for m,c in purity(t).items():P5[m]+=126*(-1)**(5-t)*comb(5,t)*2**t*c
P5={m:c for m,c in P5.items()if c};QP=product(Q2,P5)
target={m:(bound[0]*Q4.get(m,0)-QP.get(m,0),bound[1]*Q4.get(m,0))for m in set(Q4)|set(QP)};target={m:c for m,c in target.items()if c!=Z}

# Independent exact attainment in the extension t^2+A*t+B=0, w^2=u(t).
A=(F(991,1025),F(159,1025));B=(F(15,820),F(1,820))
def ta(x,y):return (add(x[0],y[0]),add(x[1],y[1]))
def ts(x,c):return (mul(x[0],c),mul(x[1],c))
def tm(x,y):
 c=mul(x[1],y[1]);return (add(mul(x[0],y[0]),neg(mul(B,c))),add(add(mul(x[0],y[1]),mul(x[1],y[0])),neg(mul(A,c))))
def tp(x,n):
 out=(ONE,Z)
 for _ in range(n):out=tm(out,x)
 return out
T=(Z,ONE);zz=ts(T,(F(75),F(-5)));uu=(Z,Z)
for i,c in enumerate([-257552,-31773836,146668933,73578600]):uu=ta(uu,ts(tp(T,i),(F(c,45433661),F(0))))
# The nested radical is exactly a root of this quadratic factor.
assert add(mul((F(991),F(159)),(F(991),F(159))),neg(scale(B,4*1025**2)))==(F(2447347),F(310013))
# Check that the cubic expression for u equals the catalogue recovery formula.
den_u=ta(ta(ts(tm(T,zz),(F(-84),F(0))),ts(T,(F(336),F(0)))),ta(tm(zz,zz),((F(-25),F(0)),Z)))
num_u=tm(T,ta(ta(ts(tm(T,zz),(F(4),F(0))),ts(T,(F(-29),F(0)))),ta(ts(zz,(F(-1),F(0))),((F(8),F(0)),Z))))
assert tm(uu,den_u)==num_u
@lru_cache(None)
def evalmon(m):
 counts=[0]*10
 for i in range(20):counts[i%10]+=m%9;m//=9
 if any(counts[i]for i in [1,2,4,5,7,8]):return (Z,)*4
 power=counts[3]+counts[9];v=tm(tp(T,counts[6]),tm(tp(zz,counts[9]),tp(uu,power//2)))
 return v+(Z,Z)if power%2==0 else (Z,Z)+v
att=[Z]*4
for m,c in target.items():att=[add(x,mul(c,y))for x,y in zip(att,evalmon(m))]
assert att==[Z]*4
# Rational interval certificate for the real radical branch and u>0.
def ia(x,y):return (x[0]+y[0],x[1]+y[1])
def im(x,y):
 vals=[a*b for a in x for b in y];return min(vals),max(vals)
def sqrt_interval(x,digits=70):
 scale10=10**digits
 return (F(isqrt(x[0].numerator*scale10**2//x[0].denominator),scale10),F(isqrt(x[1].numerator*scale10**2//x[1].denominator)+1,scale10))
sr=sqrt_interval((F(61),F(61)));inner=sqrt_interval(ia((F(2447347),)*2,im((F(310013),)*2,sr)))
tr=im(ia(ia(inner,(F(-991),)*2),im((F(-159),)*2,sr)),(F(1,2050),)*2)
assert F('-0.012846')<tr[0]<tr[1]<F('-0.012845')
ur=(F(0),F(0));power=(F(1),F(1))
for c in [-257552,-31773836,146668933,73578600]:ur=ia(ur,im((F(c,45433661),)*2,power));power=im(power,tr)
assert ur[0]>0
zr=im(im((F(5),)*2,ia((F(15),)*2,im((F(-1),)*2,sr))),tr)
dr=ia(ia(im((F(-84),)*2,im(tr,zr)),im((F(336),)*2,tr)),ia(im(zr,zr),(F(-25),)*2))
assert dr[0]*dr[1]>0

print('PASS exact attaining radical state and positive normalization.',flush=True)

def determinant(A):
 M=[row[:]for row in A];d=ONE
 for j in range(len(M)):
  pivot=next((i for i in range(j,len(M))if M[i][j]!=Z),None)
  if pivot is None:return Z
  if pivot!=j:M[j],M[pivot]=M[pivot],M[j];d=neg(d)
  v=M[j][j];d=mul(d,v)
  for i in range(j+1,len(M)):
   f=div(M[i][j],v)
   for k in range(j+1,len(M)):M[i][k]=add(M[i][k],neg(mul(f,M[j][k])))
 return d
positivity=[];prepared=[];globalden=lcm(*(v.denominator for c in target.values() for v in c))
assert len(basis['blocks'])==len(cert['blocks'])==len(cert['basisScales'])
for B,G,sc in zip(basis['blocks'],cert['blocks'],cert['basisScales']):
 assert B['kind']==G['kind'] and B['L2']==G['L2'];n=B['dimension'];H=[[read(v)for v in row]for row in G['matrix']]
 assert len(H)==n and all(len(row)==n for row in H) and all(H[a][b]==H[b][a]for a in range(n)for b in range(n))
 minors=[determinant([r[:j]for r in H[:j]])for j in range(1,n+1)];assert all(positive(v)for v in minors)
 positivity.append({'kind':B['kind'],'L2':B['L2'],'dimension':n,'leadingPrincipalMinorSigns':[1 for value in minors]})
 channels=[];den=[]
 for cov in B['covariants']:
  D=lcm(*(F(c).denominator for comp in cov for m,pair in comp for c in pair));den.append(D);polys=[]
  assert len(cov)==B['L2']+1
  for comp in cov:
   poly={}
   for mon,pair in comp:
    assert len(mon)==4 and all(0<=i<20 for i in mon)
    c=read(pair);poly[encode(mon)]=(int(c[0]*D),int(c[1]*D))
   polys.append(poly)
  channels.append(polys)
 assert len(channels)==n
 weights=[F(w)for w in B['weights']];assert weights==[F(1,comb(B['L2'],r))for r in range(B['L2']+1)]
 adjusted=[[scale(H[a][b],1/(F(sc[a])*F(sc[b])*den[a]*den[b]))for b in range(n)]for a in range(n)]
 denominator=lcm(*(v.denominator for row in adjusted for c in row for v in c))*lcm(*(w.denominator for w in weights))
 globalden=lcm(globalden,denominator)
 prepared.append((B,n,channels,weights,adjusted))
 print('PASS exact positive matrix',B['kind'],B['L2'],'dimension',n,flush=True)
# One common denominator avoids repeated reductions of huge rational sums.
result={};allkeys=set(target);denominator=globalden
for B,n,channels,weights,adjusted in prepared:
 for r,w in enumerate(weights):
  # Form G*v first, then v^dagger*(G*v). This exact distributive expansion
  # avoids recomputing identical monomial products for each matrix entry.
  for a in range(n):
   mixed={}
   for b in range(n):
    factor=scale(adjusted[a][b],w*denominator);assert all(v.denominator==1 for v in factor)
    fa,fb=map(int,factor)
    if not(fa or fb):continue
    for m,(ba,bb)in channels[b][r].items():
     key=conjugate(m);va,vb=mixed.get(key,(0,0))
     mixed[key]=(va+fa*ba+61*fb*bb,vb+fa*bb+fb*ba)
   for m,(aa,ab)in channels[a][r].items():
    for mm,(ba,bb)in mixed.items():
     key=m+mm;va,vb=result.get(key,(0,0))
     result[key]=(va+aa*ba+61*ab*bb,vb+aa*bb+ab*ba)
 print('PASS full block expansion',B['kind'],B['L2'],flush=True)
allkeys.update(result);result={m:c for m,c in result.items()if c!=(0,0)}
assert all((v*globalden).denominator==1 for c in target.values()for v in c)
targetints={m:tuple(int(v*globalden)for v in c)for m,c in target.items()}
assert result==targetints,'Full polynomial identity failed'
print('PASS full complex-polynomial identity',len(allkeys),'monomials.',flush=True)
report={'n':9,'sector':5,'boundLatex':r'\frac{130928+6832\sqrt{61}}{1461}','bound':[str(v)for v in bound],'degree':8,'fieldRadicand':61,'fullPolynomialIdentityVerified':True,'allComplexAmplitudes':True,'floatingPointUsed':False,'nonzeroMonomials':len(target),'allMonomialsIncludingCancellations':len(allkeys),'positiveDefiniteGramBlocks':len(positivity),'maxMatrixDimension':max(b['dimension']for b in basis['blocks']),'squaredForms':sum(b['dimension']*(b['L2']+1)for b in basis['blocks']),'attainmentVerified':True,'realBranchAndPositiveNormalizationVerified':True,'positivity':positivity,'basisSha256':hashlib.sha256(basis_file.read_bytes()).hexdigest(),'certificateSha256':hashlib.sha256(cert_file.read_bytes()).hexdigest()}
(BASE/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
