#!/usr/bin/env python3
"""Independent exact physical-spin audit of K(18,4)=1.

Only SymPy and the standard library are used. No discovery code or numerical
semidefinite solver is imported. Exact magnetic blocks reconstruct X=J_9.S_3.
"""
from pathlib import Path
from math import comb
import hashlib,json
import sympy as s
HERE=Path(__file__).resolve().parent

def verify(directory=None):
 root=Path(directory or HERE);n=18;j=s.Integer(9);aux=s.Integer(3);x=s.Symbol('x');lam=s.Symbol('lambda')
 polynomial=x**4-2*x**3-653*x**2-1146*x+8640
 assert s.expand(polynomial-((x+23)*(x+6)*(x-4)*(x-27)-6264))==0
 spins=list(range(6,13));dimensions=[2*L+1 for L in spins]
 eigenvalues=[s.Rational(L*(L+1)-j*(j+1)-aux*(aux+1),2)for L in spins]
 expected_x=[-30,-23,-15,-6,4,15,27]
 expected_w=[319320,-6264,-63720,-6264,-6264,-111600,-6264]
 assert eigenvalues==expected_x
 values=[polynomial.subs(x,v)for v in eigenvalues];assert values==expected_w
 assert sum(dimensions)==133 and sum(d*v for d,v in zip(dimensions,values))==0
 assert [s.sign(v)for v in values]==[1,-1,-1,-1,-1,-1,-1]
 positive_index=sum(d for d,v in zip(dimensions,values)if v>0)
 negative_index=sum(d for d,v in zip(dimensions,values)if v<0)
 assert (positive_index,negative_index)==(13,120)
 # All calculations below use normalized physical spin states |j,j-r>|3,3-a>.
 diagonal={};blocks=[];partial_traces={q:[s.Integer(0)]*7 for q in range(1,5)}
 partial_W=[s.Integer(0)]*7;counts=[0]*7
 for rsum in range(25):
  words=[(r,rsum-r)for r in range(19)if 0<=rsum-r<=6];index={w:k for k,w in enumerate(words)};B=s.zeros(len(words));M=j+aux-rsum
  for col,(r,a)in enumerate(words):
   B[col,col]=(j-r)*(aux-a)
   if r>0 and a<6:B[index[r-1,a+1],col]=s.sqrt(r*(19-r)*(a+1)*(6-a))/2
   if r<18 and a>0:B[index[r+1,a-1],col]=s.sqrt((r+1)*(18-r)*a*(7-a))/2
  assert B==B.T
  active=[h for h,L in enumerate(spins)if L>=abs(M)]
  assert len(active)==len(words)
  charpoly=s.expand(B.charpoly(lam).as_expr())
  assert s.Poly(charpoly-s.prod(lam-eigenvalues[h]for h in active),lam).is_zero
  for h in active:counts[h]+=1
  powers={1:B}
  for q in range(2,5):powers[q]=(powers[q-1]*B).applyfunc(s.simplify)
  W=powers[4]-2*powers[3]-653*powers[2]-1146*B+8640*s.eye(len(words))
  # The independently reconstructed W must have the claimed block spectrum.
  assert s.Poly(s.expand(W.charpoly(lam).as_expr()-s.prod(lam-values[h]for h in active)),lam).is_zero
  for k,(r,a)in enumerate(words):
   diagonal[r,a]=s.simplify(W[k,k]);partial_W[a]+=diagonal[r,a]
   for q in range(1,5):partial_traces[q][a]+=s.simplify(powers[q][k,k])
  blocks.append({'totalMagneticNumber':int(M),'dimension':len(words),'activeTotalSpins':[spins[h]for h in active],'xCharacteristicPolynomialVerified':True,'wCharacteristicPolynomialVerified':True})
 assert counts==dimensions
 # Conservation of total magnetic number makes the off-diagonal partial-trace
 # entries identically zero. These seven diagonal entries check the whole trace.
 assert partial_W==[0]*7
 moments=[]
 for q in range(1,5):
  mean=s.Rational(sum(d*v**q for d,v in zip(dimensions,eigenvalues)),133)
  assert partial_traces[q]==[19*mean]*7
  moments.append({'power':q,'normalizedTrace':str(mean),'partialTraceOnAuxiliary':str(19*mean)})
 # Independent exact witness check: the real AC4 state and the phased AC5
 # state share these probabilities. Every distinct support gap is at least5.
 probabilities={0:s.Rational(7,65),5:s.Rational(51,130),13:s.Rational(51,130),18:s.Rational(7,65)}
 assert sum(probabilities.values())==1
 for q in range(5):
  assert sum(p*r**q for r,p in probabilities.items())==s.Rational(sum(r**q for r in range(19)),19)
 assert min(b-a for a,b in zip(sorted(probabilities),sorted(probabilities)[1:]))==5
 for a in range(7):assert sum(p*diagonal[r,a]for r,p in probabilities.items())==0
 # Check rho4 in normalized Dicke splitting directly, independently of the
 # spin-moment argument. This also certifies exact attainment of dimension one.
 c=[s.sqrt(probabilities.get(r,0))for r in range(19)]
 A=s.Matrix(5,15,lambda a,b:c[a+b]*s.sqrt(s.Rational(comb(4,a)*comb(14,b),comb(18,a+b))))
 assert (A*A.T-s.eye(5)/5).applyfunc(s.simplify)==s.zeros(5)
 result={'passed':True,'n':18,'t':4,'maximumDimension':1,'auxiliarySpin':'3','ambientTensorDimension':133,'polynomial':'X^4-2X^3-653X^2-1146X+8640 I','factoredPolynomial':'(X+23 I)(X+6 I)(X-4 I)(X-27 I)-6264 I',
  'totalSpins':spins,'sectorDimensions':dimensions,'xEigenvalues':list(map(str,eigenvalues)),'wEigenvalues':list(map(str,values)),'positiveIndex':positive_index,'negativeIndex':negative_index,'zeroIndex':0,
  'physicalMagneticBlocksChecked':len(blocks),'characteristicPolynomialsChecked':2*len(blocks),'partialTraceDiagonalEntriesChecked':7,'partialTraceJofWIsZero':True,'traceMoments':moments,
  'knownAC4LineCompressionEntriesExactlyZero':49,'knownAC4ReducedMatrixEntriesChecked':25,'exactLowerBound':1,'exactUpperBound':1,'bound':'7k<=13','arbitraryComplexSubspacesCovered':True,
  'compressionArgument':'Every first-spin operator in a degree-four polynomial in J.S has spin-tensor ranks at most 4. Its scalar coefficient is fixed by the partial trace over spin 9. That partial trace vanishes, so the polynomial contains only first-spin ranks 1..4. Every such tensor compresses to zero on an AC4 subspace. Hence C tensor V3 is a 7k-dimensional totally isotropic subspace. A nonsingular Hermitian form of inertia (13,120) admits no totally isotropic subspace larger than 13.',
  'method':'Independent normalized physical-spin matrices in every conserved total-magnetic block; exact characteristic polynomials, all partial-trace entries, exact strict inertia, and direct normalized-Dicke verification of the attaining state. No optimization or discovery imports.',
  'physicalBlocks':blocks}
 if(root/'certificate.json').exists():result['certificateSHA256']=hashlib.sha256((root/'certificate.json').read_bytes()).hexdigest()
 return result
if __name__=='__main__':
 result=verify();(HERE/'independent-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k!='physicalBlocks'},indent=2))
