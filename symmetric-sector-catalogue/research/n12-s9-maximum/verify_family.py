"""Exact verification of a one-parameter family attaining S9=1280 (SymPy)."""
import sympy as s,json
from pathlib import Path
from math import comb
z=s.symbols('z',real=True);r=s.sqrt(462);den=27*z*z+64
A=-9*r*(81*z**4+318*z*z+64)
B=42*den+3*s.I*r*z*(27*z*z+88)
C=2*r*(27*z*z+32)+42*s.I*z*den
D=-9*den*(3*z*z+2)-72*s.I*r*z
x={0:A,2:B,4:C,6:D,8:-C,10:B,12:-A}
q=s.expand(sum(comb(12,i)*v*s.conjugate(v)for i,v in x.items()));print('q',s.factor(q),flush=True)
p=[]
for t in range(7):
 val=0
 for a in range(t+1):
  for b in range(t+1):
   rho=s.expand(sum(comb(12-t,j)*x.get(a+j,0)*s.conjugate(x.get(b+j,0))for j in range(13-t)))
   val+=comb(t,a)*comb(t,b)*s.expand(rho*s.conjugate(rho))
 p.append(s.expand(val));print('purity',t,flush=True)
p+=p[5::-1];P9=s.expand(comb(12,9)*sum((-1)**(9-t)*comb(9,t)*2**t*p[t]for t in range(10)));assert s.expand(P9-1280*q*q)==0
flip=s.expand(sum((-1)**i*comb(12,i)*x.get(i,0)*x.get(12-i,0)for i in range(13)));assert flip==0
multipoles={}
for L in [1,3,5]:
 a=s.expand(s.Rational((2*L+1)*comb(12,L),(L+1)*comb(13+L,12))*sum((-1)**(t+L)*comb(L,t)*comb(t+L,t)*p[t]for t in range(L+1)));assert a==0;multipoles[str(L)]='0'
values=[s.factor(p[2].subs(z,k)/q.subs(z,k)**2)for k in [0,1]];assert values[0]!=values[1]
report={'purity2AtParameters0And1':list(map(str,values)),'rotationInequivalentMembersVerified':True,'n':12,'sector':9,'value':'1280','parameter':'z real','binomialScaledAmplitudes':{str(i):str(s.expand(v))for i,v in x.items()},'normSquared':str(s.factor(q)),'exactSectorIdentityVerified':True,'oddMultipoles':multipoles,'spinFlipOverlap':'0','globalOptimalityEstablishedByThisCheck':False}
Path(__file__).with_name('family-verification.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS exact one-parameter family',flush=True)
