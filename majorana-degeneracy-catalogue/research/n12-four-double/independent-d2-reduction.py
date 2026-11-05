from pathlib import Path
import sympy as s,json
z,A,B,c,d,v=s.symbols('z A B c d v',real=True)
Q=z**4+(A+s.I*B)*z*z+1
P=s.Poly(s.expand(Q**2*(z**4+1+(c+s.I*d)*z*z)),z);p=[P.nth(r) for r in range(13)]
Z=s.expand(sum(q*s.conjugate(q)/s.binomial(12,r) for r,q in enumerate(p)))
R00=sum(s.binomial(9,b)*p[b]*s.conjugate(p[b])/(s.binomial(12,b)**2)for b in range(10))
R02=sum(s.binomial(9,b)*p[b]*s.conjugate(p[b+2])/(s.binomial(12,b)*s.binomial(12,b+2))for b in range(10))
E=[]
for expr in [R00-Z/4,s.re(R02),s.im(R02)]:
 expr=s.Poly(s.expand(expr),c,d);assert expr.coeff_monomial(c*d)==0 and expr.coeff_monomial(c*c)==expr.coeff_monomial(d*d)
 E.append(s.Poly(expr.coeff_monomial(c*c)*v+expr.coeff_monomial(c)*c+expr.coeff_monomial(d)*d+expr.coeff_monomial(1),A,B,c,d,v).clear_denoms()[1].primitive()[1].as_expr())
M,Y=s.linear_eq_to_matrix(E,[c,d,v]);det=s.factor(M.det());print('linear determinant',det,flush=True)
sol=M.inv()*Y;sol=[s.factor(q)for q in sol];difference=s.factor(sol[2]-sol[0]**2-sol[1]**2);print('rankone defect',difference,flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps({'variables':['A','B','c','d','v'],'equations':list(map(str,E)),'linearDeterminant':str(det),'linearSolution':list(map(str,sol)),'rankOneDefect':str(difference)},indent=2)+'\n')
