from search import *
x,a=sy.symbols('x a');f=[x**3,x-x**6/sy.Integer(3),1+3*x**5];M=sy.Matrix([[sy.diff(v,x,j).subs(x,a)for v in f]for j in range(2)]);u=[(-1)**j*M[:,[i for i in range(3)if i!=j]].det()for j in range(3)];p=sy.expand(sum(v*w for v,w in zip(u,f)));q=sy.div(p,(x-a)**2,x);assert q[1]==0
Q=sy.Poly(q[0],x);qfun=sy.lambdify(a,Q.all_coeffs(),'mpmath');pfun=sy.lambdify(a,sy.Poly(p,x).all_coeffs()[::-1],'mpmath')
def bar(t):
 roots=mp.polyroots(qfun(t),maxsteps=500);points=[xyz(t)]*2+[xyz(z)for z in roots];b=mp.matrix([sum(v[j]for v in points)/6 for j in range(3)]);return b,roots

def dist(t):b,_=bar(t);return mp.sqrt(sum(z*z for z in b))
t=mp.mpf(3520147021340202)/10**15;b,roots=bar(t);d=dist(t)
co=pfun(t);c=mp.matrix([v/((-1)**r*mp.sqrt(mp.binomial(6,r)))for r,v in enumerate(co)]);c/=mp.sqrt(sum(abs(v)**2 for v in c));c*=mp.exp(-1j*mp.arg(c[3]));P=Problem(6);cc=P.MV.T*c
R=dict(n=6,dimension=3,subspaceId=P.row['id'],distance=mp.nstr(d,80),barycenter=[mp.nstr(v,80)for v in b],basisCoordinates=[cstr(v)for v in cc],dickeCoefficients=[cstr(v)for v in c],roots=[{'root':cstr(t),'multiplicity':2}]+[{'root':cstr(r),'multiplicity':1}for r in roots],multiplicityPattern=[2,1,1,1,1],method='Stationary point on a real double-root curve; high-precision refinement of unrestricted and double-root DE searches.',stationaryResidual=mp.nstr(abs(mp.diff(dist,t)),15),repeatedRoot=mp.nstr(t,80),symbolicPolynomialBasis=list(map(str,f)),symbolicKernel=list(map(str,u)),symbolicQuotient=str(Q.as_expr()),optimalityProved=False)
R['method']='Exact rational parameter in the double-root curve, chosen near the numerical stationary point';R['exactRepeatedRoot']='3520147021340202/1000000000000000';R['exactPolynomialDivisionVerified']=True;R['exactMultiplicityTwoVerified']=bool(Q.as_expr().subs({x:sy.Rational(3520147021340202,10**15),a:sy.Rational(3520147021340202,10**15)})!=0);save('n6-rational-result.json',R);print(json.dumps(R,indent=2))
