#!/usr/bin/env python3
"""Independent exact checks of the structured creation/annihilation parent formulas."""
from pathlib import Path
import sympy as s,json,hashlib
HERE=Path(__file__).resolve().parent
simp=lambda x:s.simplify(s.expand_complex(x))

def spin(n):
 j=s.Rational(n,2);Z=s.diag(*[j-r for r in range(n+1)]);P=s.zeros(n+1)
 for r in range(1,n+1):P[r-1,r]=s.sqrt(r*(n-r+1))
 return Z,P

def maps(n,kind):
 if kind=='creation':
  A=s.zeros(n+2,n+1);B=s.zeros(n+2,n+1)
  for r in range(n+1):A[r,r]=s.sqrt(n-r+1);B[r+1,r]=s.sqrt(r+1)
 else:
  A=s.zeros(n,n+1);B=s.zeros(n,n+1)
  for r in range(n+1):
   if r<n:A[r,r]=s.sqrt(n-r)
   if r>0:B[r-1,r]=s.sqrt(r)
 return A,B

def zero(M):return all(simp(x)==0 for x in M)
def identities(n,kind):
 j=s.Rational(n,2);Z,P=spin(n);I=s.eye(n+1);A,B=maps(n,kind);Zp,Pp=spin(n+(1 if kind=='creation'else -1))
 if kind=='creation':
  checks=[A.H*A-((j+1)*I+Z),B.H*B-((j+1)*I-Z),A.H*B-P.H,
    A.H*Zp*A-((j+1)*I+Z)*(Z+I/2),B.H*Zp*B-((j+1)*I-Z)*(Z-I/2),
    A.H*Zp*B-P.H*(Z-I/2),B.H*Zp*A-P*(Z+I/2),
    A.H*Pp*A-P*((j+2)*I+Z),B.H*Pp*B-P*((j+1)*I-Z),
    A.H*Pp*B-((j+1)**2*I-Z**2),B.H*Pp*A-P**2]
 else:
  checks=[A.H*A-(j*I+Z),B.H*B-(j*I-Z),A.H*B-P,
    A.H*Zp*A-(j*I+Z)*(Z-I/2),B.H*Zp*B-(j*I-Z)*(Z+I/2),
    A.H*Zp*B-P*(Z+I/2),B.H*Zp*A-P.H*(Z-I/2),
    A.H*Pp*A-P*(j*I+Z),B.H*Pp*B-P*((j-1)*I-Z),
    A.H*Pp*B-P**2,B.H*Pp*A-(j*j*I-Z**2)]
 assert all(zero(X)for X in checks),(n,kind)
 return len(checks)

def main():
 count=0
 for n in range(2,21):
  for kind in ['creation','derivative']:count+=identities(n,kind)
 source=HERE.parent/'n18-t2/basis.json';raw=source.read_bytes();data=json.loads(raw);n=18;k=4;j=s.Integer(9)
 U=s.Matrix.hstack(*[s.Matrix([s.sympify(x)for x in b['exactCoefficients']])for b in data['basis']])
 Z,P=spin(n);c=j*(j+1)/3
 assert all(zero(X)for X in [U.H*U-s.eye(k),U.H*Z*U,U.H*P*U,U.H*Z**2*U-c*s.eye(k),U.H*P**2*U,U.H*(Z*P+P*Z)*U])
 records=[]
 for kind in ['creation','derivative']:
  A,B=maps(n,kind);beta=j+(1 if kind=='creation'else 0);E=s.Matrix.hstack(A*U,B*U)/s.sqrt(beta)
  Zp,Pp=spin(n+(1 if kind=='creation'else -1));alpha=s.Rational(n+(3 if kind=='creation'else -1),3)
  expectedZ=s.diag(*([alpha/2]*k+[-alpha/2]*k));expectedP=s.zeros(2*k)
  for a in range(k):expectedP[a,a+k]=alpha if kind=='creation'else 0;expectedP[a+k,a]=alpha if kind=='derivative'else 0
  assert zero(E.H*E-s.eye(2*k))
  assert zero(E.H*Zp*E-expectedZ)
  assert zero(E.H*Pp*E-expectedP)
  records.append({'kind':kind,'parentDimension':2*k,'physicalQubits':n+(1 if kind=='creation'else -1),'alpha':str(alpha),'orientation':'positive spinor'if kind=='creation'else 'transpose spinor; equivalent to negative spinor by a y half-turn','exactlyVerified':True})
 report={'passed':True,'operatorIdentityChecks':count,'operatorIdentitySizes':[2,20],'sourceAC2ConditionsIndependentlyChecked':True,'sourceBasisSHA256':hashlib.sha256(raw).hexdigest(),'exactN18ParentChecks':records,'proofScope':'Universal proof follows by Schwinger operator identities; these finite checks independently audit their normalization and implementation. Reproduces existing upper bounds, not a stronger one.'}
 (HERE/'conditional-parent-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
