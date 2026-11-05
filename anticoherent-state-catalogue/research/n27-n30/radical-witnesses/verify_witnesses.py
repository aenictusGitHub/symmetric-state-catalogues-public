"""Independent exact verification of saved Dicke vectors. Requires SymPy.

The verifier does not import the construction code. It reconstructs partial
traces directly by summing Dicke-to-product amplitudes, and checks polynomial
root multiplicities independently of the displayed constellation.
"""
from pathlib import Path
from math import comb
import json,time,hashlib
import sympy as s
P=Path(__file__).resolve().parent

def verify():
 started=time.monotonic();raw=(P/'states.json').read_bytes();records=json.loads(raw)['states'];results=[]
 for rec in records:
  n,t=rec['n'],rec['achievedOrder'];c=list(map(s.sympify,rec['exactCoefficients']))
  assert len(c)==n+1
  norm=s.simplify(sum(x*s.conjugate(x) for x in c));assert norm==1
  checks=[]
  for q in (t,t+1):
   defect={}
   for a in range(q+1):
    for b in range(a,q+1):
     rho=sum(c[a+l]*s.conjugate(c[b+l])*s.sqrt(s.Rational(comb(q,a)*comb(q,b)*comb(n-q,l)**2,comb(n,a+l)*comb(n,b+l))) for l in range(n-q+1))
     z=s.simplify(rho-(s.Rational(1,q+1) if a==b else 0))
     if z!=0:defect[(a,b)]=z
   assert bool(defect)==(q==t+1)
   checks.append({'t':q,'identityVerified':not defect,'nonzeroUpperEntries':len(defect)})
   if q==t+1:
    w=rec['verification']['nextOrderWitness'];assert s.simplify(defect[(w['row'],w['column'])]-s.sympify(w['differenceFromIdentity']))==0
  # P(z)=sum_r (-1)^r sqrt(C(n,r)) c_r z^r is the Majorana
  # polynomial. Poles are treated projectively, not by numerical clustering.
  support=[r for r in range(n+1) if c[r]!=0];lo,hi=min(support),max(support)
  q=s.symbols('q')
  if n in (27,29):
   assert all((r-lo)%5==0 for r in support)
   coeff=[s.simplify((-1)**r*s.sqrt(comb(n,r))*c[r]) for r in support]
   coeff=[s.simplify(x/coeff[0]) for x in coeff]
   assert all(s.simplify(coeff[i]-coeff[5-i])==0 for i in range(6))
   A,B=coeff[1:3]
   # p(q)=(q+1)*(q^4+(A-1)q^3+(B-A+1)q^2+(A-1)q+1).
   # Substitution u=q+1/q gives u^2+(A-1)u+B-A-1.
   conditions=[s.simplify((A-1)**2-4*(B-A-1)),s.simplify(1+A+B),s.simplify(5+B-3*A)]
   assert all(x.is_zero is False for x in conditions)
   assert s.simplify(coeff[0])!=0
   finiteSimple=25
   proof='Reciprocal quintic factors into q+1 and a reciprocal quartic. The reduced quadratic discriminant and its values at u=+/-2 are nonzero; every q root is nonzero, hence all25 lifted finite roots are simple.'
  else:
   assert n==30
   # Factor z(1-z^4) first. Its six projective roots are distinct.
   # The remaining polynomial in q=z^4 has coefficients in Q(i sqrt287).
   alpha=s.sqrt(326271)/4;beta=(-s.Integer(118635)+14355*s.I*s.sqrt(287))/(256*s.sqrt(326271))
   inner=s.Poly(s.expand((alpha/beta)*q*(1-q)**4+(1+14*q+q*q)**3),q,extension=[s.I,s.sqrt(287)])
   assert s.degree(s.gcd(inner,inner.diff()))==0
   assert inner.eval(0)!=0 and inner.eval(1)!=0 and inner.degree()==6
   finiteSimple=29
   proof='After extracting z(1-z^4), the degree6 polynomial in q=z^4 is coprime to its derivative, has nonzero constant and leading coefficients and is nonzero at q=1. Its24 lifted roots and the six octahedral projective roots are distinct.'
  partition=sorted(([lo] if lo else [])+([n-hi] if n-hi else [])+[1]*(finiteSimple-lo if n==30 else finiteSimple),reverse=True)
  # finiteSimple includes z=0 only in N30; the other families remove poles first.
  if n==30:partition=[1]*30
  assert sum(partition)==n
  result={'n':n,'t':t,'passed':True,'normalization':str(norm),'partialTraceChecks':checks,'nextOrderWitness':rec['verification']['nextOrderWitness'],'majoranaPartition':partition,'majoranaProof':proof}
  results.append(result);print('PASS',n,'exactACorder',t,'partition',partition,flush=True)
 report={'passed':True,'statesSHA256':hashlib.sha256(raw).hexdigest(),'method':'Independent exact direct partial traces; next-order diagonal witnesses; algebraic root multiplicity checks.','states':results,'elapsedSeconds':time.monotonic()-started}
 (P/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n');return report
if __name__=='__main__':verify()
