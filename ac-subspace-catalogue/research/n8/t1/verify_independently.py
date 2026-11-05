#!/usr/bin/env python3
"""Independent direct-3j audit of the pure16 obstruction used for N8,t1.
No 6j symbols or rational highest-weight summation routine is used here.
"""
from pathlib import Path
from functools import lru_cache
import json,hashlib,sys
import sympy as s
from sympy.physics.wigner import wigner_3j
ROOT=Path(__file__).resolve().parent
DEFAULT=ROOT/'certificate.json'

@lru_cache(None)
def tensor(L,q,m,mp):
 if -m+q+mp:return s.S.Zero
 return s.sqrt(2*L+1)*(-1)**(8-m)*wigner_3j(8,L,8,-m,q,mp)

@lru_cache(None)
def eigen(L,J):
 # Direct expectation of sum_q T_Lq tensor T_Lq^dagger in |J,J>.
 weights={m:s.sqrt(2*J+1)*wigner_3j(8,8,J,m,J-m,-J)for m in range(J-8,9)}
 assert s.simplify(sum(v*v for v in weights.values()))==1
 total=s.S.Zero
 for m,w in weights.items():
  for mp,wp in weights.items():
   q=m-mp
   if abs(q)<=L:total+=w*wp*tensor(L,q,m,mp)*tensor(L,q,J-mp,J-m)
 return s.simplify(total)

def verify(path=DEFAULT):
 raw=Path(path).read_bytes();cert=json.loads(raw);identity=cert['identity'];levels=list(range(1,8))+[9]
 for L in levels:
  for q in range(-L,L+1):
   norm=sum(tensor(L,q,m,mp)**2 for m in range(-8,9)for mp in range(-8,9))
   assert s.simplify(norm)==1,(L,q,'tensor norm')
 coeff=[s.Rational(x)for x in identity['lowMultipoleCoefficients']];channels=[]
 for J in range(0,17,2):
  lhs=s.Rational(identity['constant'])+(s.Rational(identity['p0Coefficient'])if J==0 else 0)
  lhs+=s.Rational(identity['a9Coefficient'])*eigen(9,J)
  lhs+=sum(coeff[L-1]*eigen(L,J)for L in range(1,8))
  assert s.simplify(lhs)==0,J
  channels.append({'J':J,'identityExactlyZero':True})
 assert all(s.Rational(identity[k])>0 for k in['constant','p0Coefficient','a9Coefficient'])
 # An independent simple k3 witness, without the packaged pair-module basis.
 Z=s.diag(*[4-r for r in range(9)]);Jp=s.zeros(9)
 for r in range(1,9):Jp[r-1,r]=s.sqrt(r*(9-r))
 B=s.zeros(9,3);B[0,0]=B[8,0]=s.sqrt(2)/2;B[2,1]=B[6,1]=s.sqrt(2)/2;B[4,2]=1
 assert s.simplify(B.H*B)==s.eye(3)
 assert s.simplify(B.H*Z*B)==s.zeros(3) and s.simplify(B.H*Jp*B)==s.zeros(3)
 # Intertwining identity is checked at the infinitesimal-generator level,
 # separately from the case verifier's 17 highest-weight vector checks.
 Y=s.zeros(9)
 for r in range(9):Y[r,8-r]=(-1)**r
 Jm=Jp.T;Jx=(Jp+Jm)/2;Jy=(Jp-Jm)/(2*s.I)
 for J in[Z,Jx,Jy]:assert s.simplify(J*Y+Y*J.T)==s.zeros(9)
 assert Y*Y.H==s.eye(9)
 return {'passed':True,'certificateSha256':hashlib.sha256(raw).hexdigest(),'recouplingMethod':'Direct Wigner3j tensor matrices and Clebsch-Gordan highest-weight vectors; no6j or rational-sum evaluator imported.','recouplingEntriesIndependentlyChecked':len(levels)*9,'tensorComponentsNormalizedExactly':sum(2*L+1 for L in levels),'channelChecks':channels,'preferredObstruction':identity['underSevenAC'],'spin4TimeReversalIntertwinerGeneratorsChecked':3,'separateSimpleDimension3WitnessVerified':True,'conceptualAudit':['The Casimir-map spectral inequality gives k<=j. If k=j, equality forces P=alpha I+X with X entirely in tensor rankN.','The unitary time-reversal matrix Y satisfies JY+YJ^T=0, so vectorization XY carries the adjoint rankN representation into total spinN, the fully symmetric2N-qubit sector.','For Hermitian X, the left N-qubit reduction of normalized vec(XY) is X²/Tr(X²)=I/(N+1)+4X/[N(N+2)]. No real-coefficient assumption is used.','The SU2-covariant one-qubit partial trace kills rankN because its target operator space has ranks at mostN−1. It maps I/(N+1) to I/N, giving a pure2N-qubit(N−1)-AC state.','ForN8 the exact positive pure16 identity is impossible. The verified dimension3 witness therefore attains the maximum.']}

if __name__=='__main__':
 report=verify(Path(sys.argv[1])if len(sys.argv)>1 else DEFAULT);(ROOT/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
