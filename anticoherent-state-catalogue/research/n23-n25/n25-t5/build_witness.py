"""Exact radical five-anticoherent N25 state. Requires SymPy only."""
from pathlib import Path
import json
import sympy as s
HERE=Path(__file__).resolve().parent

def reduction(c,n,t):
 M=s.Matrix(t+1,n-t+1,lambda a,b:c[a+b]*s.sqrt(s.Rational(s.binomial(t,a)*s.binomial(n-t,b),s.binomial(n,a+b))))
 return (M*M.H).applyfunc(s.simplify)

def main():
 n=25;w=(s.Integer(1001)-s.I*s.sqrt(1561327))/(16*s.sqrt(10013))
 assert s.simplify(w*s.conjugate(w))==1
 raw={0:s.sqrt(s.Rational(21,250)),25:-s.sqrt(s.Rational(21,250)),5:-s.I*s.sqrt(s.Rational(31,125)),20:s.I*s.sqrt(s.Rational(31,125)),10:-s.I*w*s.sqrt(s.Rational(21,125)),15:s.I*w*s.sqrt(s.Rational(21,125))}
 c=[s.simplify(raw.get(r,s.S.Zero)) for r in range(26)]
 assert s.simplify(sum(x*s.conjugate(x) for x in c))==1
 checks=[]
 for t in range(1,7):
  D=(reduction(c,n,t)-s.eye(t+1)/s.Integer(t+1)).applyfunc(s.simplify)
  zero=D==s.zeros(t+1);assert zero==(t<=5)
  dist=s.simplify(s.trace(D*D.H));checks.append({'t':t,'identityVerified':zero,'frobeniusDistanceSquared':str(dist),'purity':str(s.simplify(dist+s.Rational(1,t+1)))})
  if t==6:assert D[0,0]==s.Rational(13,35420)
 formula=r'\sqrt{\frac{21}{250}}\bigl(|D_{25}^{(0)}\rangle-|D_{25}^{(25)}\rangle\bigr)-i\sqrt{\frac{31}{125}}\bigl(|D_{25}^{(5)}\rangle-|D_{25}^{(20)}\rangle\bigr)-iw\sqrt{\frac{21}{125}}\bigl(|D_{25}^{(10)}\rangle-|D_{25}^{(15)}\rangle\bigr)'
 state={'id':'radical-n25-t5','n':25,'achievedOrder':5,'witnessesOrders':list(range(1,6)),'label':'Exact dihedral state · order 5','status':'proven-existence','stateLatex':formula,'parametersLatex':r'w=\frac{1001-i\sqrt{1561327}}{16\sqrt{10013}},\quad|w|=1','support':sorted(raw),'exactCoefficients':[str(x) for x in c],'coefficients':[[float(s.re(x)),float(s.im(x))] for x in c],'normalization':'sum_r |c_r|^2 = 1, verified exactly','coefficientConvention':'Normalized Dicke basis; r counts excitations.','provenance':'Rational moment matching on six C5-separated Dicke indices with an exact phase-triangle cancellation.','notes':'Exactly fifth-order AC, not sixth-order. Global maximum order requires a separate nonexistence theorem.','verification':{'method':'Exact Schmidt matrices through order six.','exactNormalization':True,'exactOrderVerified':True,'globalMaximumOrderProven':False,'reductions':checks,'nextOrderWitness':{'order':6,'row':0,'column':0,'differenceFromIdentity':'13/35420','differenceLatex':r'\frac{13}{35420}'}}}
 (HERE/'states.json').write_text(json.dumps({'states':[state]},indent=2,ensure_ascii=False)+'\n')
 print('PASS N25: rho5=I6/6; rho6[0,0]-1/7=13/35420.')
if __name__=='__main__':main()
