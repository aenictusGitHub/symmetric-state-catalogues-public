"""Exact algebraic stationarity and independent high-precision full-space checks.
Requires SymPy and mpmath. This is not a global-minimum proof.
"""
from algebraic import *
from pathlib import Path
from math import comb
import json
import mpmath as mp

u=s.symbols('u')
lo=s.Rational('0.34031561845452539230789192544421338373007474986524')
hi=s.Rational('0.34031561845452539230789192544421338373007474986525')
assert s.Poly(R,c).count_roots(lo,hi)==1
assert s.Poly(R,c).is_irreducible
print('PASS: irreducible degree14 polynomial; exactly one real root in supplied rational interval.',flush=True)
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
 return imul(x,(1/y[1],1/y[0]))
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
assert zr[0]>0 and ur[0]>0
qr=iadd(iadd((1,1),imul((120,120),ur)),iadd(imul((210,210),ipow((lo,hi),2)),imul((10,10),imul(ur,ipow(zr,2)))))
assert qr[0]>0
print('PASS: rational interval bounds certify z>0, u>0, Q>0 and nonzero recovery denominators.',flush=True)
Us=s.cancel(U.subs(z,Z));un,ud=s.fraction(Us)
assert s.gcd(ud,R)==1
# Each support equation reduces to a polynomial in c,u,z after removing b when appropriate.
Hs=[]
for f in [F[0]/b,F[1],F[2]/b]:
 di=s.Poly(s.expand(f.subs(d,z*b)),b).as_dict();H=s.expand(di.get((2,),0)*u+di.get((0,),0));Hs.append(H)
for H in Hs:
 num=s.fraction(s.cancel(H.subs(z,Z).subs(u,Us)))[0]
 assert s.rem(num,R,c)==0
print('PASS: all three exact stationary equations vanish modulo the degree14 polynomial.',flush=True)
mp.mp.dps=100
root=mp.findroot(s.lambdify(c,R,'mpmath'),(mp.mpf(str(lo.evalf(110))),mp.mpf(str(hi.evalf(110)))))
zz=s.lambdify(c,Z,'mpmath')(root);uu=s.lambdify(c,Us,'mpmath')(root)
assert uu>0
bb=mp.sqrt(uu);dd=zz*bb;QQ=1+120*uu+210*root**2+10*dd**2
coeffs=[mp.mpf('0')]*11
for i,xx in [(0,1),(3,bb),(6,root),(9,dd)]:coeffs[i]=mp.sqrt(comb(10,i))*xx/mp.sqrt(QQ)

# Independent reduced-purity matrix calculation with ambient complex derivatives.
def value_gradient(x):
 norm=mp.sqrt(sum(v*v for v in x));y=[v/norm for v in x]
 cv=[mp.mpc(y[i],y[i+11]) for i in range(11)]
 weights=[0]*6
 for t in range(7):weights[min(t,10-t)]+=comb(10,6)*(-1)**(6-t)*2**t*comb(6,t)
 value=mp.mpf('0');grad=[mp.mpc(0)]*11
 for t,w in enumerate(weights):
  if not w:continue
  fac=[[mp.sqrt(mp.mpf(comb(t,i)*comb(10-t,j))/comb(10,i+j)) for j in range(11-t)] for i in range(t+1)]
  A=mp.matrix([[cv[i+j]*fac[i][j] for j in range(11-t)] for i in range(t+1)])
  rho=A*A.H;value+=w*sum(abs(rho[i,j])**2 for i in range(t+1) for j in range(t+1))
  T=rho*A
  for i in range(t+1):
   for j in range(11-t):grad[i+j]+=4*w*T[i,j]*fac[i][j]
 gg=[v.real for v in grad]+[v.imag for v in grad];radial=sum(v*w for v,w in zip(y,gg))
 return value,[(g-v*radial)/norm for g,v in zip(gg,y)]
x=coeffs+[mp.mpf('0')]*11
value,g=value_gradient(x);gn=mp.sqrt(sum(v*v for v in g))
print('S6=',mp.nstr(value,100),flush=True);print('full complex normalized gradient norm=',mp.nstr(gn,12),flush=True)
assert gn<mp.mpf('1e-85')
result={'n':10,'k':6,'support':[0,3,6,9],'root':mp.nstr(root,100),'ratio_d_over_b':mp.nstr(zz,100),'ratio_b_squared':mp.nstr(uu,100),'coefficients':[mp.nstr(coeffs[i],100) for i in [0,3,6,9]],'sectorLength':mp.nstr(value,100),'fullComplexGradientNorm':mp.nstr(gn,30),'stationarityVerifiedExactly':True,'positiveRecoveryCertifiedByRationalIntervals':True,'globalOptimalityProved':False,'rootInterval':[str(lo),str(hi)],'rootPolynomialDegree':14}
# Numerical local Hessian only; zero directions include normalization, phase, and collective rotations.
mp.mp.dps=75;h=mp.mpf('1e-23');H=mp.matrix(22)
for j in range(22):
 xp=x.copy();xm=x.copy();xp[j]+=h;xm[j]-=h
 gp=value_gradient(xp)[1];gm=value_gradient(xm)[1]
 for i in range(22):H[i,j]=(gp[i]-gm[i])/(2*h)
H=(H+H.T)/2
evals=mp.eigsy(H,eigvals_only=True)
result['fullNormalizedHessianEigenvalues']=[mp.nstr(v,40) for v in evals]
print('Hessian eigs=',[mp.nstr(v,8) for v in evals],flush=True)
(Path(__file__).resolve().parent / 'verified.json').write_text(json.dumps(result,indent=2))
