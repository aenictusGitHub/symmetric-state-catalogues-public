"""Exact reduced-purity proof within two N=6 two-component families.
No restriction to these families is justified for a global minimum proof.
"""
from pathlib import Path
from math import comb
import sympy as s
import json
p,q,x=s.symbols('p q x',real=True)
N=6
out=[]
for support,critical in [((2,5),s.Rational(20,21)),((0,5),s.Rational(4,7))]:
 c=[s.S(0)]*(N+1);c[support[0]]=p;c[support[1]]=q
 purities=[]
 for t in range(N+1):
  A=s.Matrix(t+1,N-t+1,lambda i,j:c[i+j]*s.sqrt(s.Rational(comb(t,i)*comb(N-t,j),comb(N,i+j))))
  rho=A*A.T
  purities.append(s.expand(s.trace(rho*rho)))
 sectors=[]
 for k in range(N+1):
  F=s.expand(comb(N,k)*sum((-1)**(k-t)*comb(k,t)*2**t*purities[t] for t in range(k+1)))
  assert all(a%2==0 and b%2==0 for (a,b),coef in s.Poly(F,p,q).terms())
  normalized=s.expand(sum(coef*x**(a//2)*(1-x)**(b//2) for (a,b),coef in s.Poly(F,p,q).terms()))
  sectors.append(normalized)
 f=sectors[4];assert s.diff(f,x).subs(x,critical)==0
 value=f.subs(x,critical);assert value==s.Rational(45,7)
 weight=s.Poly(f,x).LC();assert weight>0 and 0<critical<1
 assert s.expand(f-value-weight*(x-critical)**2)==0
 achieved=[s.factor(v.subs(x,critical)) for v in sectors]
 assert sum(achieved)==64
 result={'n':6,'sector':4,'support':list(support),'firstProbability':str(critical),'secondProbability':str(1-critical),'sectorPolynomial':str(f),'squareCompletion':str(value)+' + '+str(weight)+'*(x-('+str(critical)+'))**2','allSectorPolynomials':[str(v) for v in sectors],'achievedSectorVectorIncludingS0':[str(v) for v in achieved],'globalMinimalityProvedByThisScript':False,'globalProof':'sos-certificate.json (verified separately by verify_certificate.py)','minimumWithinFamilyProved':True,'relativePhaseRemovableByCollectiveZRotation':True,'majoranaMultiplicityPattern':([2,1,1,1,1] if support==(2,5) else [1]*6),'properRotationGroup':('C3' if support==(2,5) else 'C5'),'fullPointGroup':('C3v' if support==(2,5) else 'C5v')}
 out.append(result)
 print(json.dumps(result,indent=2))
for f,g in zip(out[0]['allSectorPolynomials'],out[1]['allSectorPolynomials']):
 assert s.expand(s.sympify(f,locals={'x':x})-s.sympify(g,locals={'x':x}).subs(x,3*x/5))==0
Path(__file__).with_name('two-component-families.json').write_text(json.dumps(out,indent=2)+'\n')
