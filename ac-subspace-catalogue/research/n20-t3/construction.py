"""Explicit radical AC3 planes; coefficient recurrence used for reproduction."""
from pathlib import Path
import sympy as s,json,time
I=s.I

def ladder2(n,r,q):return s.prod((n-r-h)*(r+h+1)for h in range(q))
def n19():
 n=19;p=list(map(lambda x:s.Rational(x,500),[77,199,21,203]));R=[1,6,11,16];S=[3,8,13,18];f2=[s.sqrt(ladder2(n,r,2))for r in R];f3=[s.sqrt(ladder2(n,r,3))for r in S[:-1]]
 aa=s.simplify(f2[0]**2*p[0]*p[3]);bb=s.simplify(f2[1]**2*p[1]*p[2]);cos=s.factor((5*aa-bb)/(4*aa));w=cos+I*s.sqrt(1-cos*cos);pass
 z0=s.sqrt(aa);Z=[z0,z0*(-2+w),z0*(1-2*w),z0*w];L=s.sqrt(f3[0]**2*p[1]*p[3]);M=f3[1]*p[2];eta=-M/(2*L)+I*s.sqrt(1-M*M/(4*L*L));T=[L,M*eta,-L-M*eta]
 a=[s.sqrt(p[0])]
 for i in range(3):a.append(s.simplify(T[i]*f2[i]*a[i]/(f3[i]*Z[i])))
 b=[s.simplify(s.conjugate(Z[i]/(f2[i]*a[i])))for i in range(4)]
 return n,R,S,a,b,{'p':p,'q':p[::-1],'Z':Z,'T':T,'cos':cos,'L':L,'M':M}

def n20():
 n=20;R=list(range(0,21,4));S=list(range(3,20,4));p=[s.Rational(x,24576)for x in[3927,333,9318,5018,2483,3497]];q=[s.Rational(x,384)for x in[143,23,53,117,48]]
 U=s.Rational(89947,262144);V=s.Rational(-274443,1048576);W=s.Rational(26039,393216);Delta=s.factor(4*U*W-V*V);A=s.sqrt(U);B=(V+I*s.sqrt(Delta))/(2*A)
 Z=[(-1)**i*s.binomial(4,i)*(A+i*B)for i in range(5)];f1=[s.sqrt(ladder2(n,r,1))for r in S];f3=[s.sqrt(ladder2(n,r,3))for r in R[:-1]];Ls=[s.sqrt(s.factor(f3[i]**2*p[i]*q[i]))for i in range(5)];C=-Ls[0]+Ls[1]+Ls[2];real3=-(C*C+Ls[3]**2-Ls[4]**2)/(2*C);T3=real3+I*s.sqrt(Ls[3]**2-real3*real3);T=[-Ls[0],Ls[1],Ls[2],T3,-C-T3]
 a=[s.sqrt(p[0])]
 for i in range(5):
  expr=f3[i]*Z[i]*a[i]/(f1[i]*T[i]);a.append(s.factor(expr))
 b=[s.factor(s.conjugate(Z[i]/(f1[i]*a[i+1])))for i in range(5)]
 return n,R,S,a,b,{'p':p,'q':q,'U':U,'V':V,'W':W,'Delta':Delta,'Z':Z,'T':T,'L':Ls,'C':C}
