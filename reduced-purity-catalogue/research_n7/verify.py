import sympy as s
from sympy.physics.wigner import clebsch_gordan
import numpy as np,json
from pathlib import Path
j=s.Integer(2);T=[]
for L in range(1,5):
 block=[]
 for q in range(-L,L+1):
  block.append(s.Matrix(5,5,lambda i,k:s.sqrt(s.Rational(2*L+1,5))*clebsch_gordan(j,L,j,j-k,q,j-i)))
 T.append(block)
for block in T:
 for t in block:assert s.simplify(s.trace(t.H*t))==1
A=list(map(s.Rational,['15/14','13/14','5/7','3/7']));h=list(map(s.Rational,['-53/420','-3/28','-9/140','5/84']))
c=[0,1/s.sqrt(2),0,0,0,0,1/s.sqrt(2),0]
def marginal(c,t):
 B=s.Matrix(t+1,8-t,lambda a,b:c[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(7-t,b)/s.binomial(7,a+b)))
 return s.simplify(B*B.H)
r3=marginal(c,3);r4=marginal(c,4);phi=s.diag(0,0,1,0,0);delta=r4-(s.eye(5)-phi)/4
assert r3==s.diag(s.Rational(2,7),s.Rational(3,14),s.Rational(3,14),s.Rational(2,7))
assert s.trace(r3*r3)==s.Rational(25,98)
for L,block in enumerate(T):
 project=lambda z:sum((s.trace(t.H*z)*t for t in block),s.zeros(5))
 assert s.simplify(project(delta)+h[L]*project(phi))==s.zeros(5)
# Verify the exact coefficient reduction used in the proof.
x,y=s.symbols('x y');aa=[x/10,(y-2)/14,(x-5*y+30)/40,(-35*x+15*y+54)/280]
f=list(map(s.Rational,['9/8','7/8','1/2','0']));tau=s.Rational(4,7);mu=s.Rational(5,98)
k=[s.factor(tau*z/16-(mu-tau*z/2)**2/(4*(1+tau*(z-1)))) for z in f]
assert all(s.simplify(1+tau*(z-1)-a)==0 for z,a in zip(f,A))
assert all(s.simplify((mu-tau*z/2)/(2*a)-b)==0 for z,a,b in zip(f,A,h))
assert s.simplify(sum(b*z for b,z in zip(k,aa))-s.Rational(1,196)-s.Rational(13,411600)*(91*x-15*y+90))==0
# Independent numerical marginal/projection checks on arbitrary complex states.
tensors=[[np.array(t,dtype=complex) for t in b] for b in T];rng=np.random.default_rng(73125);worst=0
B={t:np.array([[[float(s.sqrt(s.binomial(t,a)*s.binomial(7-t,b)/s.binomial(7,a+b))) if a+b==k else 0 for b in range(8-t)] for a in range(t+1)] for k in range(8)]) for t in [3,4]}
for _ in range(300):
 z=rng.normal(size=8)+1j*rng.normal(size=8);z/=np.linalg.norm(z)
 rr={}
 for t in [3,4]:
  b=np.einsum('k,kab->ab',z,B[t]);rr[t]=b@b.conj().T
 w,v=np.linalg.eigh(rr[4]);v=v[:,0];ph=np.outer(v,v.conj());d=rr[4]-(np.eye(5)-ph)/4
 powers=[];sq=0
 for L,block in enumerate(tensors):
  cp=np.array([np.trace(t.conj().T@ph) for t in block]);cd=np.array([np.trace(t.conj().T@d) for t in block]);powers.append(np.vdot(cp,cp).real);sq+=float(A[L])*np.vdot(cd+float(h[L])*cp,cd+float(h[L])*cp).real
 xx=10*powers[0];yy=2+14*powers[1];rhs=sq+(1183*xx+195*(6-yy))/411600
 lhs=np.trace(rr[3]@rr[3]).real-25/98;worst=max(worst,abs(lhs-rhs));assert abs(lhs-rhs)<1e-11
report={'exact_attainment':True,'exact_coefficients_and_remainder':True,'arbitrary_complex_states_checked':300,'maximum_identity_residual':worst,'bound':'25/98'}
Path(__file__).with_name('verification.json').write_text(json.dumps(report,indent=2));print(report)
