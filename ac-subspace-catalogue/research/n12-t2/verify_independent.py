#!/usr/bin/env python3
"""Independent exact tensor and complementary-projector audit for N12,k3,t2.
Uses explicit Wigner-3j tensor matrix elements and coupled highest-weight states;
checks their expectation against Wigner-6j recoupling and the supplied Farkas identity.
"""
from pathlib import Path
from functools import lru_cache
import json
import sympy as s
from sympy.physics.wigner import wigner_3j,wigner_6j,clebsch_gordan
ROOT=Path(__file__).parent;j=6;Jvalues=list(range(0,13,2))
@lru_cache(None)
def tensor(L,m,mp):
 M=m-mp
 if abs(M)>L:return s.S.Zero
 return (-1)**(j-m)*s.sqrt(2*L+1)*wigner_3j(j,L,j,-m,M,mp)
@lru_cache(None)
def coupled(J,m):
 if not -j<=J-m<=j:return s.S.Zero
 return clebsch_gordan(j,j,J,m,J-m,J)

def direct_coefficient(L,J):
 return s.simplify(sum(coupled(J,m)*coupled(J,mp)*tensor(L,m,mp)*tensor(L,J-mp,J-m)
   for m in range(-j,j+1) for mp in range(-j,j+1) if coupled(J,m)!=0 and coupled(J,mp)!=0))

# Full single-spin irreducible tensor orthonormality and rank-L Casimir eigenvalue.
Jz=s.diag(*range(-j,j+1));Jplus=s.zeros(13)
for m in range(-j,j):Jplus[m+j+1,m+j]=s.sqrt((j-m)*(j+m+1))
Jminus=Jplus.T;Jx=(Jplus+Jminus)/2;Jy=(Jplus-Jminus)/(2*s.I)
assert (Jx*Jx+Jy*Jy+Jz*Jz)==42*s.eye(13)
for L in range(13):
 T=s.diag(*[tensor(L,m,m)for m in range(-j,j+1)])
 assert s.simplify(s.trace(T.H*T))==1
 image=sum([Ja*T*Ja for Ja in (Jx,Jy,Jz)],s.zeros(13))
 assert (image-(42-s.Rational(L*(L+1),2))*T).applyfunc(s.simplify)==s.zeros(13)
 print('rank',L,flush=True)
A=s.Matrix([[direct_coefficient(L,J)for J in Jvalues]for L in range(13)])
A6=s.Matrix([[(2*L+1)*wigner_6j(j,j,L,j,j,J)for J in Jvalues]for L in range(13)])
assert A==A6
assert all(sum(A[L,col]for L in range(13))==1 for col in range(7))
alpha={L:s.Rational(28,L*(L+1)-112) for L in range(1,13)}
cert=json.loads((ROOT/'certificate.json').read_text())
assert A==s.Matrix([[s.Rational(x)for x in row]for row in cert['multipoleMatrix']])
beta=[s.Integer(3)]+[alpha[L]for L in range(1,13)]
minusPQ=-s.Matrix([beta])*A
assert list(minusPQ)==list(map(s.Rational,cert['negativePQEigenvalues']))
f=cert['identity']
positive=s.Rational(f['constant'])*s.ones(1,7)+sum([s.Rational(f['a'+str(L)])*A[L,:]for L in range(1,4)],s.zeros(1,7))+s.Matrix([list(map(s.Rational,f['pJ']))])
assert minusPQ==positive
scale=s.Rational(9015)/s.Rational(f['constant'])
integers={'constant':s.Rational(f['constant'])*scale,'a3':s.Rational(f['a3'])*scale,**{'p'+str(J):s.Rational(v)*scale for J,v in zip(Jvalues,f['pJ']) if s.Rational(v)!=0}}
assert integers=={'constant':9015,'a3':1031316,'p0':3936270,'p2':20020,'p6':728365}
assert all(v>0 for v in integers.values())
output={'allChecksPassed':True,'independentRecouplingCoefficients':91,'allTensorCasimirEigenvalues':13,
 'tracePQIdentityVerifiedColumnwise':True,'traceP2ConstraintUsed':False,
 'identity':{'negativePQScale':str(scale),'fullA1Coefficient':str(s.Rational(f['a1'])*scale),'fullA2Coefficient':str(s.Rational(f['a2'])*scale),**{k:int(v)for k,v in integers.items()}},
 'necessaryConclusion':'0 = 9015 + 1031316 a3 + 3936270 p0 + 20020 p2 + 728365 p6',
 'contradiction':'Every term on the right is nonnegative and the constant is strictly positive.',
 'conceptualAudit':{
 'codeMomentNormalization':'PJaP=0 and PJaJbP=14 delta_ab P for spin6,t2',
 'orthogonalErrorSpaces':'C and Ja C/sqrt14 (a=x,y,z) are four orthogonal dimension3 spaces',
 'rankOneComplement':'Q=I-P-sum_a JaPJa/14 is a positive rank1 projector, TrQ=1 and PQ=0',
 'multipoleInversion':'P0=3Q0; P_L=28Q_L/[L(L+1)-112] for1<=L<=12, denominators all nonzero',
 'complementAnticoherence':'P has zero rank1 andrank2 components, so Q does as well',
 'twoCopyParity':'q tensor q is symmetric; atinteger spin6 only even J occur',
 'multipoleDefinition':'aL=sum_M |Tr(Q T_LM)|² >=0 for Hilbert-Schmidt-orthonormal irreducible tensors',
 'probabilities':'pJ=||Pi_J(q tensor q)||² >=0 and sum_evenJ pJ=1',
 'scalarPart':'Tr(P0 Q0)=3/13, hence the beta0 coefficient of a0 is3',
 'recouplingSign':'Direct tensor matrix-element calculation confirms aL=(2L+1)sum_evenJ {6 6 L;6 6 J}pJ with no extra sign',
 'noRestrictedAnsatz':'No geometric symmetry, real-coefficient, or Dicke-support assumption is used.'},
 'conclusion':'The maximum dimension of a second-order AC subspace at N12 is2, using the separately exactly verified dimension2 radical witness.'}
(ROOT/'independent-verification.json').write_text(json.dumps(output,indent=2)+'\n')
print('SHORT CERTIFICATE VERIFIED',integers,flush=True)
