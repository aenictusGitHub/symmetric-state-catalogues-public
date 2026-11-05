import sympy as s,json
from pathlib import Path
R=Path(__file__).resolve().parent
n=11;k=6
a,b,c,d=s.symbols('a b c d',real=True)
xs={0:a,3:b,6:c,9:d}
Q=a*a+165*b*b+462*c*c+55*d*d
w=[0]*6
for r in range(k+1):w[min(r,n-r)]+=s.binomial(n,k)*(-1)**(k-r)*s.binomial(k,r)*2**r
F=0
for r,wr in enumerate(w):
 if not wr:continue
 # rho_{i,j}=sqrt(C(r,i)C(r,j))*sum_t C(n-r,t) x_{i+t}x_{j+t}.
 purity=0
 for i in range(r+1):
  for j in range(r+1):
   inner=sum(s.binomial(n-r,t)*xs.get(i+t,0)*xs.get(j+t,0) for t in range(n-r+1))
   purity+=s.binomial(r,i)*s.binomial(r,j)*inner**2
 F+=wr*purity
F=s.expand(F)
print('Q=',Q)
print('F=',F)
print('gradienteqs projective')
G=[]
for v in (a,b,c):
 g=s.factor((Q*s.diff(F,v)-2*F*s.diff(Q,v)).subs(d,1))
 G.append(g)
 print(v,g)
(R/'polynomials.json').write_text(json.dumps({'variables':['a','b','c','d'],'normalization':str(Q),'sectorPolynomial':str(F),'projectiveEquations':[str(x) for x in G]},indent=2)+'\n')
if __name__=='__main__':
 import mpmath as mp
 mp.mp.dps=100
 amp=[mp.mpf('-.4863546496685766'),mp.mpf('.2609151148419338'),mp.mpf('.5760723489633521'),mp.mpf('.6029287738635656')]
 initial=[amp[i]/mp.sqrt([1,165,462][i])/(amp[3]/mp.sqrt(55)) for i in range(3)]
 fn=s.lambdify((a,b,c),G,'mpmath')
 sol=mp.findroot(fn,initial,tol=mp.mpf('1e-95'),maxsteps=100)
 qq=s.lambdify((a,b,c),Q.subs(d,1),'mpmath')(*sol)
 coeff=[sol[0]/mp.sqrt(qq),mp.sqrt(165)*sol[1]/mp.sqrt(qq),mp.sqrt(462)*sol[2]/mp.sqrt(qq),mp.sqrt(55)/mp.sqrt(qq)]
 vv=s.lambdify((a,b,c),F.subs(d,1),'mpmath')(*sol)/qq**2
 print('ratios',*[mp.nstr(x,100) for x in sol],sep='\n')
 print('coefficients',*[mp.nstr(x,100) for x in coeff],sep='\n')
 print('value',mp.nstr(vv,100))
 (R/'high-precision.json').write_text(json.dumps({'ratios':[mp.nstr(x,100) for x in sol],'coefficients':[mp.nstr(x,100) for x in coeff],'value':mp.nstr(vv,100)},indent=2)+'\n')
