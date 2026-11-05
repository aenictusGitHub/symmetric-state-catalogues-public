"""Exact radical pure AC states at N27,N29,N30. Requires SymPy only."""
from pathlib import Path
import json
import sympy as s
P=Path(__file__).resolve().parent

def reduction(c,n,t):
 M=s.Matrix(t+1,n-t+1,lambda a,b:c[a+b]*s.sqrt(s.Rational(s.binomial(t,a)*s.binomial(n-t,b),s.binomial(n,a+b))))
 return (M*M.H).applyfunc(s.simplify)

def main():
 sources=[]
 for n,start,nums,den,x,y,z in [(27,1,[327,494,429],s.Integer(50),52,3433,6137),(29,2,[3598,2081,3696],25*s.sqrt(30),3094,2170247,11743083)]:
  w=(s.Integer(x)-s.I*s.sqrt(y))/s.sqrt(z);assert s.simplify(w*s.conjugate(w))==1
  a,b,c=[s.sqrt(q)/den for q in nums]
  raw={start:a,start+5:-s.I*b,start+10:-s.I*w*c,n-start-10:s.I*w*c,n-start-5:s.I*b,n-start:-a}
  formula=(r'\sqrt{\frac{'+str(nums[0])+'}{'+s.latex(den**2)+r'}}\bigl(|D_{'+str(n)+r'}^{('+str(start)+r')}\rangle-|D_{'+str(n)+r'}^{('+str(n-start)+r')}\rangle\bigr)'
    r'-i\sqrt{\frac{'+str(nums[1])+'}{'+s.latex(den**2)+r'}}\bigl(|D_{'+str(n)+r'}^{('+str(start+5)+r')}\rangle-|D_{'+str(n)+r'}^{('+str(n-start-5)+r')}\rangle\bigr)'
    r'-iw\sqrt{\frac{'+str(nums[2])+'}{'+s.latex(den**2)+r'}}\bigl(|D_{'+str(n)+r'}^{('+str(start+10)+r')}\rangle-|D_{'+str(n)+r'}^{('+str(n-start-10)+r')}\rangle\bigr)')
  sources.append((n,5,raw,formula,r'w=\frac{'+str(x)+r'-i\sqrt{'+str(y)+r'}}{\sqrt{'+str(z)+r'}}',{'family':'C5 phase-triangle','wRealNumerator':x,'wImaginaryRadicand':y,'wDenominatorRadicand':z}))
 n=30;z=s.symbols('z');F6=z*(1-z**4);F8=1+14*z**4+z**8;alpha=s.sqrt(326271)/4;beta=(-s.Integer(118635)+14355*s.I*s.sqrt(287))/(256*s.sqrt(326271));poly=s.Poly(alpha*F6**5+beta*F6*F8**3,z);raw={r:s.simplify(poly.nth(r)/s.sqrt(s.binomial(n,r))) for r in range(n+1) if poly.nth(r)!=0}
 formula=r'\frac{\beta}{\sqrt{30}}\bigl(|D_{30}^{(1)}\rangle-|D_{30}^{(29)}\rangle\bigr)+\frac{\alpha+41\beta}{\sqrt{142506}}\bigl(|D_{30}^{(5)}\rangle-|D_{30}^{(25)}\rangle\bigr)+\frac{-5\alpha+549\beta}{\sqrt{14307150}}\bigl(|D_{30}^{(9)}\rangle-|D_{30}^{(21)}\rangle\bigr)+\frac{10\alpha+2237\beta}{\sqrt{119759850}}\bigl(|D_{30}^{(13)}\rangle-|D_{30}^{(17)}\rangle\bigr)'
 sources.append((30,7,raw,formula,r'\alpha=\frac{\sqrt{326271}}4,\quad\beta=\frac{-118635+14355i\sqrt{287}}{256\sqrt{326271}}',{'family':'Octahedral binary forms','polynomial':'alpha*(z*(1-z**4))**5+beta*z*(1-z**4)*(1+14*z**4+z**8)**3','alpha':str(alpha),'beta':str(beta)}))
 records=[]
 for n,t,raw,formula,params,extra in sources:
  c=[s.simplify(raw.get(r,s.S.Zero)) for r in range(n+1)];assert s.simplify(sum(x*s.conjugate(x) for x in c))==1
  checks=[];witness=None
  for q in range(1,t+2):
   D=(reduction(c,n,q)-s.eye(q+1)/s.Integer(q+1)).applyfunc(s.simplify);zero=D==s.zeros(q+1);assert zero==(q<=t)
   dist=s.simplify(s.trace(D*D.H));checks.append({'t':q,'identityVerified':zero,'frobeniusDistanceSquared':str(dist),'purity':str(s.simplify(dist+s.Rational(1,q+1)))})
   if q==t+1:witness={'order':q,'row':0,'column':0,'differenceFromIdentity':str(D[0,0]),'differenceLatex':s.latex(D[0,0])};assert D[0,0]!=0
  records.append({'id':('octahedral-' if n==30 else 'dihedral-')+f'n{n}-t{t}','n':n,'achievedOrder':t,'witnessesOrders':list(range(1,t+1)),'label':f'Exact '+('octahedral' if n==30 else 'dihedral')+f' state · order {t}','status':'proven-existence','stateLatex':formula,'parametersLatex':params,'support':sorted(raw),'exactCoefficients':[str(x) for x in c],'coefficients':[[float(s.re(x)),float(s.im(x))] for x in c],'normalization':'sum_r |c_r|^2 = 1, verified exactly','coefficientConvention':'Normalized Dicke basis; r counts excitations.','provenance':'Independent analytical construction by '+extra['family']+'.','construction':extra,'notes':f'Exactly order {t}, not order {t+1}. Global maximum order requires a separate exclusion theorem.','verification':{'method':'Exact Schmidt matrices through the next order.','exactNormalization':True,'exactOrderVerified':True,'globalMaximumOrderProven':False,'reductions':checks,'nextOrderWitness':witness}})
  print('EXACT',n,t,'nextorderdefect',witness['differenceFromIdentity'],flush=True)
 (P/'states.json').write_text(json.dumps({'states':records},indent=2,ensure_ascii=False)+'\n')
if __name__=='__main__':main()
