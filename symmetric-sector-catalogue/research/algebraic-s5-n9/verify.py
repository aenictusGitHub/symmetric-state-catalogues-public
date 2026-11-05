"""Exact algebraic stationarity and independent high-precision full-space checks.
Requires SymPy and mpmath. This is not a global-maximum proof.
"""
from radical import *
N=9
K=5
DIM=N+1
from pathlib import Path
from math import comb
import json
import mpmath as mp

u=s.symbols('u')
lo=s.Rational('-0.01284555902806194669499550925759084638053248768653')
hi=s.Rational('-0.01284555902806194669499550925759084638053248768652')
assert s.Poly(R,c).count_roots(lo,hi)==1
assert s.Poly(R,c).is_irreducible
print('PASS: irreducible degree4 polynomial; exactly one real root in supplied rational interval.',flush=True)
zn,zd=s.fraction(Z)
assert s.gcd(zd,R)==1
# Rational interval arithmetic certifies the recovered branch is real and positive.
def iadd(x,y):return (x[0]+y[0],x[1]+y[1])
def imul(x,y):
 vals=[a*b for a in x for b in y];return (min(vals),max(vals))
def ipow(x,p):
 ans=(s.Integer(1),s.Integer(1))
 for _ in range(p):ans=imul(ans,x)
 return ans
def idiv(x,y):
 assert y[0]*y[1]>0
 return imul(x,(s.Rational(1)/y[1],s.Rational(1)/y[0]))
def ipoly(expr,gens,ranges):
 ans=(s.Integer(0),s.Integer(0))
 for powers,coef in s.Poly(expr,*gens).terms():
  term=(coef,coef)
  for power,ran in zip(powers,ranges):term=imul(term,ipow(ran,power))
  ans=iadd(ans,term)
 return ans
zr=idiv(ipoly(zn,[c],[(lo,hi)]),ipoly(zd,[c],[(lo,hi)]))
un0,ud0=s.fraction(U)
ur=idiv(ipoly(un0,[c,z],[(lo,hi),zr]),ipoly(ud0,[c,z],[(lo,hi),zr]))
assert zr[1]<0 and ur[0]>0
qr=iadd(iadd((1,1),imul((84,84),ur)),iadd(imul((84,84),ipow((lo,hi),2)),imul(ur,ipow(zr,2))))
assert qr[0]>0
print('PASS: rational interval bounds certify z<0, u>0, Q>0 and nonzero recovery denominators.',flush=True)
# Prove that the displayed nested radical chooses the same isolated root.
from math import isqrt
def isqrt_interval(x,digits=90):
 scale=10**digits
 low,high=map(s.Rational,x)
 assert low>=0
 ll=isqrt(int(s.floor(low*scale**2)))
 hh=isqrt(int(s.floor(high*scale**2)))+1
 return (s.Rational(ll,scale),s.Rational(hh,scale))
sr=isqrt_interval((61,61))
ir=isqrt_interval(iadd((2447347,2447347),imul((310013,310013),sr)))
tr=idiv(iadd(iadd((-991,-991),imul((-159,-159),sr)),ir),(2050,2050))
assert lo<tr[0]<=tr[1]<hi
assert s.simplify(QUADRATIC.subs(c,T))==0
assert s.rem(R,QUADRATIC,c,extension=S)==0
assert s.rem(Z-ZSHORT,QUADRATIC,c,extension=S)==0
print('PASS: compact radical t is the isolated quartic root and z=5(15-sqrt(61))t exactly.',flush=True)
Us=s.cancel(U.subs(z,Z));un,ud=s.fraction(Us)
assert s.gcd(ud,R)==1
# Each support equation reduces to a polynomial in c,u,z after removing b when appropriate.
Hs=[]
for f in [F[0]/b,F[1],F[2]/b]:
 di=s.Poly(s.expand(f.subs(d,z*b)),b).as_dict();H=s.expand(di.get((2,),0)*u+di.get((0,),0));Hs.append(H)
for H in Hs:
 num=s.fraction(s.cancel(H.subs(z,Z).subs(u,Us)))[0]
 assert s.rem(num,R,c)==0
print('PASS: all three exact stationary equations vanish modulo the degree4 polynomial.',flush=True)
# Exact attained value, using Euler's identity at a stationary point.
V=126*(1+84*u-4*c*c-u*z*z-320*u*c+80*u*c*z)/(1+84*u+84*c*c+u*z*z)
VV=s.cancel(V.subs(z,Z).subs(u,Us));vn,vd=s.fraction(VV)
assert s.rem(vn-NU*vd,QUADRATIC,c,extension=S)==0
print('PASS: exact S5=(130928+6832*sqrt(61))/1461.',flush=True)
mp.mp.dps=100
root=mp.findroot(s.lambdify(c,R,'mpmath'),(mp.mpf(str(lo.evalf(110))),mp.mpf(str(hi.evalf(110)))))
zz=s.lambdify(c,Z,'mpmath')(root);uu=s.lambdify(c,Us,'mpmath')(root)
assert uu>0
bb=mp.sqrt(uu);dd=zz*bb;QQ=1+84*uu+84*root**2+dd**2
coeffs=[mp.mpf('0')]*DIM
for i,xx in [(0,1),(3,bb),(6,root),(9,dd)]:coeffs[i]=mp.sqrt(comb(N,i))*xx/mp.sqrt(QQ)

# Independent reduced-purity matrix calculation with ambient complex derivatives.
def value_gradient(x):
 norm=mp.sqrt(sum(v*v for v in x));y=[v/norm for v in x]
 cv=[mp.mpc(y[i],y[i+DIM]) for i in range(DIM)]
 weights=[0]*(N//2+1)
 for t in range(K+1):weights[min(t,N-t)]+=comb(N,K)*(-1)**(K-t)*2**t*comb(K,t)
 value=mp.mpf('0');grad=[mp.mpc(0)]*DIM
 for t,w in enumerate(weights):
  if not w:continue
  fac=[[mp.sqrt(mp.mpf(comb(t,i)*comb(N-t,j))/comb(N,i+j)) for j in range(DIM-t)] for i in range(t+1)]
  A=mp.matrix([[cv[i+j]*fac[i][j] for j in range(DIM-t)] for i in range(t+1)])
  rho=A*A.H;value+=w*sum(abs(rho[i,j])**2 for i in range(t+1) for j in range(t+1))
  T=rho*A
  for i in range(t+1):
   for j in range(DIM-t):grad[i+j]+=4*w*T[i,j]*fac[i][j]
 gg=[v.real for v in grad]+[v.imag for v in grad];radial=sum(v*w for v,w in zip(y,gg))
 return value,[(g-v*radial)/norm for g,v in zip(gg,y)]
x=coeffs+[mp.mpf('0')]*DIM
value,g=value_gradient(x);gn=mp.sqrt(sum(v*v for v in g))
print('S5=',mp.nstr(value,100),flush=True);print('full complex normalized gradient norm=',mp.nstr(gn,12),flush=True)
assert gn<mp.mpf('1e-85')
assert abs(value-(130928+6832*mp.sqrt(61))/1461)<mp.mpf('1e-90')
result={'n':N,'k':K,'support':[0,3,6,9],'root':mp.nstr(root,100),'ratio_d_over_b':mp.nstr(zz,100),'ratio_b_squared':mp.nstr(uu,100),'coefficients':[mp.nstr(coeffs[i],100) for i in [0,3,6,9]],'sectorLength':mp.nstr(value,100),'fullComplexGradientNorm':mp.nstr(gn,30),'stationarityVerifiedExactly':True,'positiveRecoveryCertifiedByRationalIntervals':True,'globalOptimalityProved':False,'rootInterval':[str(lo),str(hi)],'rootPolynomialDegree':4,'radicalBranchVerifiedExactly':True,'exactValue':str(NU),'exactValueLatex':r'\frac{130928+6832\sqrt{61}}{1461}'}
# Numerical local Hessian only; zero directions include normalization, phase, and collective rotations.
mp.mp.dps=75;h=mp.mpf('1e-23');H=mp.matrix(2*DIM)
for j in range(2*DIM):
 xp=x.copy();xm=x.copy();xp[j]+=h;xm[j]-=h
 gp=value_gradient(xp)[1];gm=value_gradient(xm)[1]
 for i in range(2*DIM):H[i,j]=(gp[i]-gm[i])/(2*h)
H=(H+H.T)/2
evals=mp.eigsy(H,eigvals_only=True)
result['fullNormalizedHessianEigenvalues']=[mp.nstr(v,40) for v in evals]
print('Hessian eigs=',[mp.nstr(v,8) for v in evals],flush=True)
# Global optimality is supplied by a separate exact SOS proof, not by the
# local stationarity/Hessian calculation in this script.
proof_directory=Path(__file__).resolve().parent.parent/'n9-s5-maximum'
proof_report=json.loads((proof_directory/'independent-verification.json').read_text())
import hashlib
assert proof_report['certificateSha256']==hashlib.sha256((proof_directory/'sos-certificate.json').read_bytes()).hexdigest()
assert proof_report['fullPolynomialIdentityVerified'] and proof_report['attainmentVerified']
result.update(globalOptimalityProved=True,globalOptimalityVerifiedByThisScript=False,globalProofUrl='../n9-s5-maximum/index.html')
(Path(__file__).resolve().parent / 'verified.json').write_text(json.dumps(result,indent=2))
