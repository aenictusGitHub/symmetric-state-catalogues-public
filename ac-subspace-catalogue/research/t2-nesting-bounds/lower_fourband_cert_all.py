import sympy as s,json
from pathlib import Path
q,i,u,z=s.symbols('q i u z',nonnegative=True)

def objects(rho):
 n=8*q+rho;j=n/2;h=2*q+(2 if rho>=6 else 0);offset=(2+rho//2) if rho<6 else (rho-2)//2
 x=(j-offset-2*i)/h;v=j*(j+1)/(3*h*h)
 W=[v/4+x*x/4-s.Rational(35,36)*x+s.Rational(5,6),-v/4-x*x/4+s.Rational(5,12)*x+s.Rational(1,2),-v/4-x*x/4+s.Rational(13,12)*x-s.Rational(1,2),v/4+x*x/4-s.Rational(19,36)*x+s.Rational(1,6)]
 if rho==5:W=[w+s.Rational(1,2)/q*(x-s.Rational(3,2))*a for w,a in zip(W,[-s.Rational(1,3),1,-1,s.Rational(1,3)])]
 assert s.factor(sum(W)-1)==0
 assert s.factor(sum((x-l)*W[l]for l in range(4)))==0
 assert s.factor(sum((x-l)**2*W[l]for l in range(4))-v)==0
 den=576*(q+(1 if rho>=6 else 0))**2
 W=[s.cancel(den*w)for w in W];assert all(s.denom(w)==1 for w in W)
 E=[]
 for l in range(4):
  rr=offset+2*i+l*h
  E.append(W[l]*W[l].subs(i,i+1)*(rr+1)*(rr+2)*(n-rr)*(n-rr-1))
 return W,E

def poly(expr,subs):return s.Poly(s.expand(expr.subs(subs,simultaneous=True)),u,z)
def check(rho):
 W,E=objects(rho);rows=[];edges=[]
 rowcharts=[{i:a,q:z+5}for a in range(3)]+[{i:u+3,q:u+z+5}]
 edgecharts=[{i:a,q:z+5}for a in range(2)]+[{i:u+2,q:u+z+5}]
 for sub in rowcharts:
  Wp=[poly(w,sub)for w in W]
  assert all(c>=0 for p in Wp for c in p.coeffs()),('weights',rho,sub)
  assert all(p.coeff_monomial(1)>0 for p in Wp)
  rows.append({'substitution':str(sub),'terms':[len(w.terms())for w in Wp],'positiveConstants':[str(w.coeff_monomial(1))for w in Wp]})
 for sub in edgecharts:
  Ep=[poly(e,sub)for e in E]
  assert all(c>=0 for p in Ep for c in p.coeffs()),('E',rho,sub)
  F=sum(Ep)**2-2*sum(e*e for e in Ep)
  neg=[(m,c)for m,c in F.terms()if c<0]
  A=s.Poly(sum(-c*u**m[0]*z**m[1]for m,c in neg),u,z)
  P=64*Ep[0]*Ep[1]*Ep[2]*Ep[3]-A*A
  assert all(c>=0 for c in P.coeffs()),('product',rho,sub,[(m,c)for m,c in P.terms()if c<0][:2])
  assert P.coeff_monomial(1)>0
  edges.append({'substitution':str(sub),'A':str(A.as_expr()),'FplusAterms':len((F+A).terms()),'productTerms':len(P.terms()),'positiveConstants':[str((F+A).coeff_monomial(1)),str(P.coeff_monomial(1))]})
 print('rho',rho,'PASS',flush=True)
 return {'rho':rho,'rowCharts':rows,'edgeCharts':edges}
if __name__=='__main__':
 rows=[check(rho)for rho in range(8)]
 Path(__file__).with_name('lower-fourband-strengthened-certificate.json').write_text(json.dumps({'passed':True,'theorem':'K(N,2)>=floor(N/8)-1 for every integer N>=40','arithmetic':'Exact integer polynomial coefficient certificates','records':rows},indent=2)+'\n')
