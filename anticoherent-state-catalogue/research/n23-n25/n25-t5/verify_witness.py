"""Independent direct partial traces from saved radical coefficients. SymPy only."""
from pathlib import Path
import json,hashlib
import sympy as s

def main():
 here=Path(__file__).resolve().parent;path=here/'states.json';reports=[]
 for st in json.loads(path.read_text())['states']:
  n,t=st['n'],st['achievedOrder'];c=list(map(s.sympify,st['exactCoefficients']));assert(n,t)==(25,5)
  assert s.simplify(sum(x*s.conjugate(x) for x in c))==1
  support=[r for r,x in enumerate(c) if x!=0];checks=[]
  for q in range(1,7):
   R=s.zeros(q+1)
   for a in range(q+1):
    for b in range(a,q+1):
     val=0
     for r in support:
      ell=r-a
      if 0<=ell<=n-q and 0<=b+ell<=n:
       val+=c[r]*s.conjugate(c[b+ell])*s.binomial(n-q,ell)*s.sqrt(s.Rational(s.binomial(q,a)*s.binomial(q,b),s.binomial(n,r)*s.binomial(n,b+ell)))
     R[a,b]=s.simplify(val);R[b,a]=s.conjugate(R[a,b])
   D=(R-s.eye(q+1)/s.Integer(q+1)).applyfunc(s.simplify)
   zero=D==s.zeros(q+1);assert zero==(q<=5)
   dist=s.simplify(s.trace(D*D.H));assert s.simplify(dist-s.sympify(st['verification']['reductions'][q-1]['frobeniusDistanceSquared']))==0
   checks.append({'t':q,'identityVerified':zero,'frobeniusDistanceSquared':str(dist)})
   if q==6:assert D[0,0]==s.Rational(13,35420)
  z,q=s.symbols('z q');R=s.sqrt(39215);S=s.sqrt(1561327)
  A=2*s.I*R;B=-R*(S+s.I*1001)/124
  P=1+A*q+B*q*q+B*q**3+A*q**4+q**5
  actual=sum((-1)**r*s.sqrt(s.binomial(n,r))*x*z**r/c[0] for r,x in enumerate(c))
  assert s.simplify(actual-P.subs(q,z**5))==0
  # P=(q+1) q^2 [u^2+(A-1)u+B-A-1], u=q+1/q.
  # Distinct u roots, neither +/-2, give four distinct quartic roots.
  disc=s.expand((A-1)**2-4*(B-A-1))
  assert s.simplify(s.im(disc))!=0
  assert s.simplify(s.im(1+A+B))!=0
  assert s.simplify(s.im(5+B-3*A))!=0
  reports.append({'id':st['id'],'n':n,'exactOrder':t,'normalizationVerified':True,'majoranaStarsAllDistinct':True,'reductions':checks})
 (here/'independent-verification.json').write_text(json.dumps({'passed':True,'statesSha256':hashlib.sha256(path.read_bytes()).hexdigest(),'method':'Independent exact direct partial-trace sums and exact reciprocal-polynomial discriminants; no construction imports.','records':reports},indent=2)+'\n')
 print('Independent PASS: exactly5AC at N25; all25Majorana stars distinct.')
if __name__=='__main__':main()
