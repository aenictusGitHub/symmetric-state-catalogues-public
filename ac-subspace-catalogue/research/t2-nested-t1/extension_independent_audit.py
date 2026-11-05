#!/usr/bin/env python3
"""Independent full-matrix audit of the N18 containing-AC1 obstruction."""
from pathlib import Path
import json,hashlib
import sympy as s
ROOT=Path(__file__).resolve().parent
simp=lambda z:s.simplify(s.expand_complex(z))
clean=lambda A:A.applyfunc(simp)
def main():
 path=ROOT.parent/'n18-t2/basis.json';data=json.loads(path.read_text())
 B=s.Matrix.hstack(*[s.Matrix([s.sympify(v)for v in u['exactCoefficients']])for u in data['basis']])
 n=18;j=9;c=30;eye=s.eye(n+1);Z=s.diag(*[j-r for r in range(n+1)])
 U=s.zeros(n+1)
 for r in range(1,n+1):U[r-1,r]=s.sqrt(r*(n-r+1))
 X=(U+U.H)/2;Y=(U-U.H)/(2*s.I)
 assert clean(B.H*B)==s.eye(4)
 for A in (X,Y,Z):assert clean(B.H*A*B)==s.zeros(4)
 for a,A in enumerate((X,Y,Z)):
  for b,C in enumerate((X,Y,Z)):
   assert clean(B.H*A*C*B)==c*int(a==b)*s.eye(4)
 # Orthogonal error-frame construction independent of componentwise recurrence.
 F=B.row_join(X*B/s.sqrt(c)).row_join(Y*B/s.sqrt(c)).row_join(Z*B/s.sqrt(c))
 assert clean(F.H*F)==s.eye(16)
 R=clean(eye-F*F.H)
 Q=s.Matrix.hstack(*[R[:,r]/s.sqrt(R[r,r])for r in (1,0,2)])
 assert clean(Q.H*Q)==s.eye(3)
 assert clean(Q*Q.H)==R
 assert clean(F.H*Q)==s.zeros(16,3)
 D=s.diag(1,1,-s.I);T=Q*D
 Splus=s.Matrix([[0,s.sqrt(2),0],[0,0,s.sqrt(2)],[0,0,0]])
 Sx=(Splus+Splus.H)/2;Sy=(Splus-Splus.H)/(2*s.I);Sz=s.diag(1,0,-1)
 alpha=s.Rational(9,2)
 assert clean(T.H*X*T)==alpha*Sx
 assert clean(T.H*Y*T)==-alpha*Sy
 assert clean(T.H*Z*T)==alpha*Sz
 # An explicit one-dimensional AC1 addition, including all cross moments.
 v=T[:,1];Parent=B.row_join(v)
 assert clean(Parent.H*Parent)==s.eye(5)
 for A in (X,Y,Z):assert clean(Parent.H*A*Parent)==s.zeros(5)
 report={'passed':True,'basisSHA256':hashlib.sha256(path.read_bytes()).hexdigest(),
   'inputFirstAndAllOrderedSecondMomentsRechecked':True,'errorFrameGramIdentity16':True,
   'complementRank':3,'complementProjectorFromFullErrorFrame':True,
   'complementBasisColumns':[1,0,2],'phaseCorrection':['1','1','-I'],
   'compressedSpinTriple':'(9/2)*(Sx,-Sy,Sz)',
   'orientationCaveat':'Reflection of y axis; zero-compression equivalent to spin1, not an oriented su2 representation.',
   'largestAdditionDimension':1,'exactFiveDimensionalAC1ParentVerified':True,
   'largestAC1ParentDimension':5,'maximumUnrestrictedAC1Dimension':7,
   'noTwoDimensionalAdditionProof':'A 2D totally isotropic subspace for Sz in C3 must contain its kernel |0>. Its other vector is a|+1>+b|-1>; zero <0|Sx|v> and <0|Sy|v> force a=b=0. The kernel vector itself gives a 1D addition.',
   'scope':'This exact N18 dimension4 AC2 code, including any collective rotation of it.'}
 (ROOT/'extension-independent-audit.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({'passed':True,'maximumContainingAC1Dimension':5}))
if __name__=='__main__':main()
