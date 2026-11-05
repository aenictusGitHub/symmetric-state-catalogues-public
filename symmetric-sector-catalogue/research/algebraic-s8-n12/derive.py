import sympy as s,json,mpmath as mp
from pathlib import Path
R=Path(__file__).parent;a,b,c,d=s.symbols('a b c d',real=True);A,B,T=s.symbols('A B T',real=True);xs={1:a,4:b,7:c,10:d};n=12;k=8
Q=12*a*a+495*b*b+792*c*c+66*d*d;w=[0]*7
for t in range(k+1):w[min(t,n-t)]+=s.binomial(n,k)*(-1)**(k-t)*s.binomial(k,t)*2**t
F=0
for t,wt in enumerate(w):
 if not wt:continue
 purity=0
 for i in range(t+1):
  for j in range(t+1):
   inner=sum(s.binomial(n-t,r)*xs.get(i+r,0)*xs.get(j+r,0)for r in range(n-t+1))
   purity+=s.binomial(t,i)*s.binomial(t,j)*inner**2
 F+=wt*purity
F=s.expand(F);Fd=s.diff(F,d).subs(d,1);G=[66*s.diff(F,v).subs(d,1)-weight*v*Fd for v,weight in [(a,12),(b,495),(c,792)]]
G=[s.Poly(g,a,b,c).primitive()[1].as_expr()for g in G];print('F',F,flush=True);print('G',G,flush=True)
# a=A*c, b=B, c²=T. Divide odd c factors before replacing.
reduced=[]
for g in G:
 h=s.Poly(s.expand(g.subs(a,A*c).subs(b,B)),c);power=min(m[0]for m,co in h.terms());assert power in [0,1]
 h=s.Poly(s.cancel(h.as_expr()/c**power),c);assert all(m[0]%2==0 for m,co in h.terms());reduced.append(s.expand(sum(co*T**(m[0]//2)for m,co in h.terms())))
print('reduced',reduced,flush=True);assert s.degree(reduced[2],T)==1
V=s.cancel(s.solve(reduced[2],T)[0]);num,den=s.fraction(V)
h=[s.Poly(s.cancel(g.subs(T,V)*den),A,B).primitive()[1].as_expr()for g in reduced[:2]]
res=s.factor(s.resultant(h[0],h[1],A));print('resultant',res,flush=True)
mp.mp.dps=100;amps=list(map(mp.mpf,['.2698136900462214','.5511240512973504','.6979890897422282','.36915319766467286']));initial=[amps[i]/mp.sqrt([12,495,792][i])/(amps[3]/mp.sqrt(66))for i in range(3)];root=mp.findroot(s.lambdify((a,b,c),G,'mpmath'),initial,tol=mp.mpf('1e-95'))
An=root[0]/root[2];Bn=root[1];print('B',mp.nstr(Bn,90),flush=True)
factors=s.factor_list(res)[1];poly=min([f for f,m in factors],key=lambda f:abs(s.lambdify(B,f,'mpmath')(Bn)));print('selected degree',s.degree(poly,B),flush=True)
subs=s.subresultants(h[0],h[1],A);linear=next(f for f in reversed(subs)if s.degree(f,A)==1);ar=s.cancel(-linear.subs(A,0)/s.diff(linear,A));an,ad=s.fraction(ar);an=s.rem(an,poly,B);ad=s.rem(ad,poly,B);ar=s.cancel(an/ad)
qq=s.lambdify((a,b,c),Q.subs(d,1),'mpmath')(*root);value=s.lambdify((a,b,c),F.subs(d,1),'mpmath')(*root)/qq**2;coeff=[root[i]*mp.sqrt([12,495,792][i])/mp.sqrt(qq)for i in range(3)]+[mp.sqrt(66)/mp.sqrt(qq)]
D={'n':12,'k':8,'support':[1,4,7,10],'normalization':str(Q),'sectorPolynomial':str(F),'stationaryEquations':[str(g)for g in G],'reducedEquations':[str(g)for g in reduced],'Bpolynomial':str(poly),'Anumerator':str(an),'Adenominator':str(ad),'Tnumerator':str(num),'Tdenominator':str(den),'h1':str(h[0]),'h2':str(h[1]),'Bdecimal':mp.nstr(Bn,100),'Binterval':[str(int(mp.floor(Bn*10**60)))+'/10^60',str(int(mp.floor(Bn*10**60))+1)+'/10^60'],'coefficients100':[mp.nstr(v,100)for v in coeff],'valueDecimal100':mp.nstr(value,100)}
(R/'exact-algebraic.json').write_text(json.dumps(D,indent=2)+'\n');print('value',mp.nstr(value,100),flush=True)
