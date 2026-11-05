import sympy as s,json,time
from pathlib import Path
R=Path(__file__).resolve().parent
A,B,T=s.symbols('A B T',real=True)
g1=T*A**3-16*T*A**2*B-19*A*B**2+352*T*A*B-110*T*A-31*A-880*B**2+880*B
g2=T*(3*A**2*B-48*A*B**2-32*A*B+16*A+1056*B**2-178*B-352)-241*B**3+91*B
g3=T*(A**2-84*A*B+1848*B+917)-10*A*B**2+10*A*B-706*B**2-440*B-39
den=A**2-84*A*B+1848*B+917
num=10*A*B**2-10*A*B+706*B**2+440*B+39
h1=s.factor(g1.subs(T,num/den)*den)
h2=s.factor(g2.subs(T,num/den)*den)
print('h1',h1,flush=True);print('h2',h2,flush=True)
res=s.factor(s.resultant(h1,h2,A))
print('resultant B',res,flush=True)
(R/'reduced-polynomials.json').write_text(json.dumps({'variables':['A','B','T'],'A':'a/c','B':'b','T':'c^2','equations':[str(g1),str(g2),str(g3)],'Tnumerator':str(num),'Tdenominator':str(den),'h1':str(h1),'h2':str(h2),'resultantB':str(res)},indent=2))
