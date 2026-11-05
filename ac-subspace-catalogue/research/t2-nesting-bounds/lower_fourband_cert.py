import sympy as s,json
from pathlib import Path
q,i,z=s.symbols('q i z',nonnegative=True)

def cert(rho):
 n=8*q+rho;j=n/2;h=2*q+(2 if rho>=6 else 0);offset=(2+rho//2) if rho<6 else (rho-2)//2
 x=(j-offset-2*i)/h;v=j*(j+1)/(3*h*h)
 W=[v/4+x*x/4-s.Rational(35,36)*x+s.Rational(5,6),-v/4-x*x/4+s.Rational(5,12)*x+s.Rational(1,2),-v/4-x*x/4+s.Rational(13,12)*x-s.Rational(1,2),v/4+x*x/4-s.Rational(19,36)*x+s.Rational(1,6)]
 if rho==5:
  W=[w+s.Rational(1,2)/q*(x-s.Rational(3,2))*a for w,a in zip(W,[-s.Rational(1,3),1,-1,s.Rational(1,3)])]
 den=576*(q+(1 if rho>=6 else 0))**2
 W=[s.cancel(den*w)for w in W]
 assert all(s.denom(w)==1 for w in W)
 # First try all q>=3, extending down to q2 for rows.
 WP=[s.Poly(s.expand(w.subs(q,i+z+2)),i,z)for w in W]
 print('rho',rho,'Wnegative',[sum(bool(c<0) for c in p.coeffs())for p in WP],flush=True)
 E=[]
 for l in range(4):
  rr=offset+2*i+l*h
  e=W[l]*W[l].subs(i,i+1)*(rr+1)*(rr+2)*(n-rr)*(n-rr-1)
  E.append(s.Poly(s.expand(e.subs(q,i+z+3)),i,z))
 F=sum(E)**2-2*sum(e*e for e in E)
 neg=[(m,c)for m,c in F.terms()if c<0]
 print('Fnegative',len(neg),'terms',len(F.terms()),flush=True)
 A=s.Poly(sum(-c*i**m[0]*z**m[1]for m,c in neg),i,z)
 P=64*E[0]*E[1]*E[2]*E[3]-A*A
 bad=[(m,c)for m,c in P.terms()if c<0]
 print('productnegative',len(bad),'min?',bad[:2],flush=True)
 return {'rho':rho,'negativeWeightCoefficients':[sum(bool(c<0) for c in p.coeffs())for p in WP],'negativeFterms':len(neg),'negativeProductTerms':len(bad),'A':str(A.as_expr()),'works':all(c>=0 for p in WP for c in p.coeffs())and not bad}
if __name__=='__main__':
 rows=[]
 for rho in range(8):rows.append(cert(rho))
 Path(__file__).with_name('lower-fourband-certificate-attempt.json').write_text(json.dumps(rows,indent=2)+'\n')
