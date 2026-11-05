import json,sympy as s,mpmath as mp
from pathlib import Path
R=Path(__file__).resolve().parent
A,B,T=s.symbols('A B T',real=True)
raw=json.loads((R/'exact-algebraic.json').read_text())
loc={'A':A,'B':B,'T':T}
p=s.Poly(s.sympify(raw['Bpolynomial'],locals=loc),B)
num=s.sympify(raw['Anumerator'],locals=loc);den=s.sympify(raw['Adenominator'],locals=loc)
numred=s.rem(num,p.as_expr(),B);denred=s.rem(den,p.as_expr(),B)
ratio=s.cancel(numred/denred);numred,denred=s.fraction(ratio)
raw['AformulaReduced']=str(ratio);raw['AnumeratorReduced']=str(numred);raw['AdenominatorReduced']=str(denred)
lo=s.Rational(2498461150300332175000977146562933123120,10**40)
hi=s.Rational(2498461150300332175000977146562933123121,10**40)
count=p.count_roots(lo,hi)
assert count==1
assert p.is_irreducible
print('B polynomial irreducible degree',p.degree(),'roots in interval',count,flush=True)
mod71=s.factor_list(p,modulus=71)[1]
assert sorted(f.degree() for f,m in mod71)==[2,9]
assert all(m==1 for f,m in mod71)
print('Squarefree mod-71 factor degrees 2 and 9 certify the S_11 Galois group together with irreducibility.',flush=True)
def rem(f):return s.rem(f,p.as_expr(),B)
abar=rem(num*s.invert(den,p.as_expr(),B))
tnum=s.sympify(raw['Tnumerator'],locals=loc);tden=s.sympify(raw['Tdenominator'],locals=loc)
tn=rem(tnum.subs(A,abar));td=rem(tden.subs(A,abar));tbar=rem(tn*s.invert(td,p.as_expr(),B))
graw=json.loads((R/'reduced-polynomials.json').read_text())
res=[]
for text in graw['equations']:
 g=s.Poly(s.sympify(text,locals=loc),A,B,T)
 rr=0
 for (ia,ib,it),coef in g.terms():
  rr=rem(rr+coef*rem(s.pow_mod(abar,ia,p.as_expr(),B) if hasattr(s,'pow_mod') else rem(abar**ia))*B**ib*rem(tbar**it))
 assert rr==0
 res.append(str(rr))
print('All three exact stationary equations reduce to zero modulo P(B).',flush=True)

# Verify the full normalized complex gradient at 100-digit precision.
mp.mp.dps=110
bb=mp.findroot(s.lambdify(B,p.as_expr(),'mpmath'),mp.mpf('.2498461150300332175'))
aa=s.lambdify(B,ratio,'mpmath')(bb)
tt=s.lambdify((A,B),tnum/tden,'mpmath')(aa,bb)
Q=aa**2*tt+165*bb**2+462*tt+55
cc=[mp.mpf(0)]*12
for i,v in zip((0,3,6,9),(aa*mp.sqrt(tt),mp.sqrt(165)*bb,mp.sqrt(462*tt),mp.sqrt(55))):cc[i]=v/mp.sqrt(Q)
n=11;k=6;weights=[0]*6
for r in range(k+1):weights[min(r,n-r)]+=int(s.binomial(n,k))*((-1)**(k-r))*int(s.binomial(k,r))*2**r
grad=[mp.mpf(0)]*12;value=mp.mpf(0)
for r,w in enumerate(weights):
 if not w:continue
 M=mp.matrix(r+1,n-r+1)
 for i in range(r+1):
  for j in range(n-r+1):M[i,j]=cc[i+j]*mp.sqrt(mp.binomial(r,i)*mp.binomial(n-r,j)/mp.binomial(n,i+j))
 rho=M*M.T;rm=rho*M
 value+=w*sum(rho[i,j]**2 for i in range(r+1) for j in range(r+1))
 for i in range(r+1):
  for j in range(n-r+1):grad[i+j]+=4*w*rm[i,j]*mp.sqrt(mp.binomial(r,i)*mp.binomial(n-r,j)/mp.binomial(n,i+j))
residual=max(abs(grad[i]-4*value*cc[i]) for i in range(12))
assert residual<mp.mpf('1e-95')
print('Full normalized real gradient residual',mp.nstr(residual,10),'imaginary gradient exactly zero for this real state.',flush=True)
report={'rootCount':int(count),'rootInterval':[str(lo),str(hi)],'polynomialIrreducible':True,'degree':11,'stationaryRemainders':res,'fullComplexStationarityExplanation':'For real amplitudes, all imaginary-direction derivatives vanish by conjugation invariance. C3 residue support forces off-support real derivatives to vanish exactly. On-support derivatives vanish by the exact polynomial remainders.','fullRealGradientResidual':mp.nstr(residual,10),'coefficients100digits':{str(i):mp.nstr(cc[i],100) for i in (0,3,6,9)},'ratios100digits':{'A':mp.nstr(aa,100),'B':mp.nstr(bb,100),'T':mp.nstr(tt,100)},'normalization100digits':mp.nstr(Q,100),'sectorValue100digits':mp.nstr(value,100),'globalOptimalityProven':False}
report['galoisGroup']='S_11'
report['mod71Factorization']=[{'factor':str(f.as_expr()),'degree':f.degree(),'multiplicity':int(m)} for f,m in mod71]
value_data=json.loads((R/'value-polynomial.json').read_text())
V=s.symbols('V')
vp=s.Poly(s.sympify(value_data['polynomial'],locals={'V':V}),V)
vlo=s.Rational('113.1090435663466909908720436191490453102178')
vhi=s.Rational('113.1090435663466909908720436191490453102179')
assert vp.count_roots(vlo,vhi)==1
assert vp.is_irreducible
assert mp.mpf(str(vlo.p))/int(vlo.q)<value<mp.mpf(str(vhi.p))/int(vhi.q)
report['sectorValueIsolatingInterval']=[str(vlo),str(vhi)]
report['sectorValuePolynomialRootCount']=1
print('Sector value is the unique root of its degree-11 polynomial in the stated narrow rational interval.',flush=True)
(R/'exact-algebraic.json').write_text(json.dumps(raw,indent=2)+'\n')
(R/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
print('Reduced A numerator',numred,flush=True);print('Reduced A denominator',denred,flush=True)
