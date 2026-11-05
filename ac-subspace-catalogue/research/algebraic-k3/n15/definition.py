import sympy as s
from math import comb
X=s.symbols('a b_re b_im c_re c_im u v_re v_im w_re w_im x_re x_im y_re y_im',real=True)
a,bR,bI,cR,cI,u,vR,vI,wR,wI,xR,xI,yR,yI=X
A={0:a,3:(bR+s.I*bI)/16,6:(cR+s.I*cI)/64,9:(cR+s.I*cI)/64,12:(bR+s.I*bI)/16,15:a}
B={1:u/4,4:(vR+s.I*vI)/32,7:(wR+s.I*wI)/64,10:(xR+s.I*xI)/64,13:(yR+s.I*yI)/8}
C={15-r:z for r,z in B.items()}
def normmoment(V,p):return s.expand(sum(comb(15,r)*z*s.conjugate(z)*s.Rational(15-2*r,2)**p for r,z in V.items()))
def plus(V,W,k,z=False):
 return s.expand(sum(comb(15,r)*s.prod(range(16-r-k,16-r))*v*s.conjugate(W[r+k])*(14-2*r if z else 1) for r,v in V.items() if r+k in W))
E=[normmoment(A,0)-1,normmoment(A,2)-s.Rational(85,4),normmoment(B,0)-1,normmoment(B,1),normmoment(B,2)-s.Rational(85,4)]
for e in [plus(A,B,1),plus(A,B,1,True),plus(B,A,2)]:E += [s.re(e).expand(),s.im(e).expand()]
E += [s.re(plus(B,C,1)).expand(),s.im(plus(B,C,1,True)).expand(),s.re(plus(C,B,2)).expand()]
E=[s.Poly(e,*X) for e in E]
# Rational scaling improves numerical condition without changing zero set.
E=[s.expand(p.as_expr()/max(abs(co) for co in p.coeffs())) for p in E]
