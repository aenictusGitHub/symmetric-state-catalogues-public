import sympy as s,json
from pathlib import Path
from sympy.physics.wigner import wigner_6j
ROOT=Path(__file__).resolve().parent
N=13;D=N+1;j=s.Rational(N,2);k=3;qrank=D-4*k;erank=3*k
I=s.eye(D)
F=s.Matrix([[(-1)**(N+J)*(2*L+1)*wigner_6j(j,j,J,j,j,L)for J in range(D)]for L in range(D)])
T=s.Matrix([[(-1)**(N+L+K)*(2*L+1)*wigner_6j(j,j,L,j,j,K)for K in range(D)]for L in range(D)])
gamma=s.diag(*([s.Rational(qrank,k)]+[s.Rational(2*L*(L+1),65)-4 for L in range(1,D)]))
delta=s.diag(*([s.Rational(erank,k)]+[-1-gamma[L,L]for L in range(1,D)]));coefs=[I,gamma,delta];ranks=[k,qrank,erank]
E=I[0:3,:].col_join(s.ones(1,D)).col_join(s.ones(1,D)*gamma).col_join(s.ones(1,D)*gamma**2).col_join(T[1:3,:]);b=s.Matrix([s.Rational(k*k,D),0,0,k,0,qrank,0,0]);A=s.zeros(0,D);names=[]
for aa,C in enumerate(coefs):
 for bb,B in enumerate(coefs[aa:],aa):
  A=A.col_join(T*C*B).col_join(F.inv()*C*B);names += [('B',aa,bb,L)for L in range(D)]+[('p',aa,bb,L)for L in range(D)]
  if aa==bb:A=A.col_join((T-I/ranks[aa])*C*C);names += [('Cauchy',aa,aa,L)for L in range(D)]
# Rank-three complement to the nine orthogonal Schmidt partners.
M=N-2;DH=M+1;jH=s.Rational(M,2);h_rank=DH-3*k
FH=s.Matrix([[(-1)**(M+J)*(2*L+1)*wigner_6j(jH,jH,J,jH,jH,L)for J in range(DH)]for L in range(DH)])
TH=s.Matrix([[(-1)**(M+L+K)*(2*L+1)*wigner_6j(jH,jH,L,jH,jH,K)for K in range(DH)]for L in range(DH)])
r2=[s.Rational((N+L+1)*(N-L)*(N+L)*(N-L-1),N*N*(N-1)*(N-1))for L in range(DH)]
h=[s.Rational(h_rank,k)*s.sqrt(s.Rational(D,DH))]+[-3*s.sqrt(r2[L])for L in range(1,DH)]
HH=s.zeros(DH,D);HM=s.zeros(DH,D)
for L in range(DH):HH[L,L]=h[L]**2;HM[L,L]=h[L]
E=E.col_join(s.ones(1,DH)*HH);b=b.col_join(s.Matrix([h_rank]))
HI=s.eye(DH);HR=s.diag(s.Rational(DH-h_rank,h_rank),*([-1]*(DH-1)))
for aa,C in enumerate([HI,HR]):
 for bb,B in enumerate([HI,HR][aa:],aa):
  A=A.col_join(TH*C*B*HH).col_join(FH.inv()*C*B*HH);names += [('BH',aa,bb,L)for L in range(DH)]+[('pH',aa,bb,L)for L in range(DH)]
  if aa==bb:A=A.col_join((TH-HI/[h_rank,DH-h_rank][aa])*C*C*HH);names += [('CauchyH',aa,aa,L)for L in range(DH)]
# Every mixed-spin tensor-product probability, with and without time reversal.
GM=s.Matrix([[(-1)**(j+jH+J)*(2*J+1)*wigner_6j(j,jH,J,jH,j,L)for L in range(DH)]for J in range(1,N)])
flip=s.diag(*[(-1)**L for L in range(DH)])
for aa,C in enumerate(coefs):
 for bb,B in enumerate([HI,HR]):
  for ff,S in enumerate([HI,flip]):
   mat=(GM*S*B*HM*C).applyfunc(s.simplify);assert all(x.is_Rational for x in mat)
   A=A.col_join(mat);names += [('mixedPH',aa,bb,ff,J)for J in range(1,N)]

point=s.Matrix([s.Rational(v)for v in json.loads((ROOT/'relaxation-report.json').read_text())['point']])
assert E*point==b and all(v>=0 for v in A*point) and all(v>=0 for v in point)
report={'passed':True,'scope':'Exact feasibility of this necessary scalar relaxation only; no physical subspace is constructed.','pointConvention':'a_L=||P_L||_HS^2 for hypothetical rank-three P','n':13,'hypotheticalDimension':3,'inequalitiesChecked':A.rows,'equalitiesChecked':E.rows,'equalityRank':E.rank(),'arithmetic':'Exact rational arithmetic, with Wigner6j radicals simplified before evaluation.'}
(ROOT/'scalar-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
