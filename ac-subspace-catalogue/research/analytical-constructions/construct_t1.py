#!/usr/bin/env python3
"""Exact first-order anticoherent subspaces from separated reflection pairs.
All output amplitudes use normalized Dicke states, r=j-m.
"""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parent

def simp(x):return s.simplify(s.expand_complex(x))

def coeff_record(v):
 return {'exactCoefficients':[str(simp(c))for c in v],
         'coefficients':[[float(s.re(c)),float(s.im(c))]for c in v],
         'support':[r for r,c in enumerate(v)if c!=0]}

def split(v,t):
 n=len(v)-1
 return s.Matrix(t+1,n-t+1,lambda r,l:v[r+l]*s.sqrt(s.binomial(t,r)*s.binomial(n-t,l)/s.binomial(n,r+l)))

def validate(n,basis):
 V=s.Matrix.hstack(*[s.Matrix(v)for v in basis]);k=len(basis)
 assert (V.H*V).applyfunc(simp)==s.eye(k)
 mats=[split(v,1)for v in basis];checks=[]
 for i,A in enumerate(mats):
  for j,B in enumerate(mats):
   assert (A*B.H-s.eye(2)*int(i==j)/2).applyfunc(simp)==s.zeros(2),(n,i,j)
   checks.append({'i':i,'j':j,'reducedCrossIdentityVerified':True})
 mats2=[split(v,2)for v in basis];nextfailure=None
 for i,A in enumerate(mats2):
  for j,B in enumerate(mats2):
   residual=(A*B.H-s.eye(3)*int(i==j)/3).applyfunc(simp)
   if residual!=s.zeros(3)and nextfailure is None:
    nextfailure={'i':i,'j':j,'squaredFrobeniusDefect':str(simp(sum(x*s.conjugate(x)for x in residual)))}
 assert nextfailure is not None
 return {'orthonormal':True,'allCrossReductionsVerified':True,'checks':checks,'scalarEntriesChecked':4*k*k,'nextOrderFailure':nextfailure}

def reflection_pairs(n):
 j=s.Rational(n,2);basis=[]
 def vector(parts):
  row=[s.S.Zero]*(n+1)
  for m,value in parts.items():row[int(j-m)]=value
  return row
 if n%2==0:
  residue=int(j)%3
  central=residue!=2
  xs=list(range(2 if central else 1,int(j),3))
  if central:basis.append(vector({s.S.Zero:s.S.One}))
 else:
  assert n%6==5
  central=False
  xs=[s.Rational(3,2)+3*l for l in range((n+1)//6)]
 for x in xs:
  x=s.Rational(x);a=s.sqrt((x+1)/(2*x+1));b=s.sqrt(x/(2*x+1))
  basis.append(vector({x:a,-x-1:b}))
  basis.append(vector({x+1:b,-x:-a}))
 assert len(basis)==n//3+1
 return {'id':f'reflection-pairs-n{n}-k{len(basis)}-t1','n':n,'dimension':len(basis),'order':1,
  'label':'Separated reflection-pair construction','basis':[coeff_record(v)for v in basis],
  'construction':{'type':'reflection-pairs','magneticPairLowerValues':list(map(str,xs)),'includeMZero':central,
   'basisConvention':'Normalized Dicke basis r=j-m; the magnetic spin basis |m> is |D_n^(j-m)>.',
   'moduleFormula':'v_x=sqrt((x+1)/(2x+1))|x>+sqrt(x/(2x+1))|-x-1>; w_x=sqrt(x/(2x+1))|x+1>-sqrt((x+1)/(2x+1))|-x>',
   'argument':'Both vectors have zero Jz expectation. Their only possible internal J+ matrix element is a sum of two equal-magnitude opposite-sign terms; spin-ladder coefficients agree by reflection. Magnetic pairs are separated by one unused level, so different modules have zero compressed Jz and J±. If included, m=0 is separated from every module by the unused m=±1 levels.'},
  'verification':validate(n,basis),'provenance':'Independent analytic construction from reflected adjacent magnetic levels.'}

def n7():
 n=7;q=(-8+s.I*s.sqrt(41))/s.sqrt(105);basis=[]
 for parts in [{0:q/2,3:-1/s.sqrt(3),6:s.sqrt(s.Rational(5,12))},
               {1:s.sqrt(s.Rational(5,12)),4:-1/s.sqrt(3),7:q/2},
               {2:1/s.sqrt(2),5:1/s.sqrt(2)}]:
  basis.append([simp(parts.get(r,0))for r in range(n+1)])
 return {'id':'cyclic-radical-n7-k3-t1','n':7,'dimension':3,'order':1,'label':'Exact three-dimensional cyclic construction',
  'basis':[coeff_record(v)for v in basis],
  'construction':{'type':'cyclic-C3','phase':str(q),'unitPhaseIdentity':'(64+41)/105=1',
    'supports':[[0,3,6],[1,4,7],[2,5]],
    'argument':'Disjoint residue classes give orthonormality and zero off-diagonal Jz. Each probability distribution has mean r=7/2. Two cross J+ sums cancel pairwise; the third is sqrt(105)*(q+conjugate(q))/12+4/3=0.'},
  'verification':validate(n,basis),'provenance':'Independent exact radical reconstruction; all cross reductions independently verified here.'}

def central_odd(n):
 j=s.Rational(n,2);X=j*(j+1)
 assert n%6 in (1,3) and n>=7
 def vector(parts):
  row=[s.S.Zero]*(n+1)
  for m,value in parts.items():row[int(j-m)]=simp(value)
  return row
 if n%6==1:
  a=8*X/(9*(4*X-7));b=7*(4*X-15)/(18*(4*X-7));c=7*(4*X-3)/(18*(4*X-7))
  cosine=-s.sqrt(X+s.Rational(1,4))*b/(2*s.sqrt((X-s.Rational(35,4))*a*c))
  q=simp(cosine+s.I*s.sqrt(1-cosine**2))
  parts=[{s.Rational(7,2):s.sqrt(a)*q,s.Rational(1,2):-s.sqrt(b),-s.Rational(5,2):s.sqrt(c)},
         {s.Rational(5,2):s.sqrt(c),-s.Rational(1,2):-s.sqrt(b),-s.Rational(7,2):s.sqrt(a)*q},
         {s.Rational(3,2):1/s.sqrt(2),-s.Rational(3,2):1/s.sqrt(2)}]
  first_x=s.Rational(11,2);kind='central-three'
 else:
  a=16*X/(52*X-135);b=27*(4*X-35)/(8*(52*X-135));c=45*(4*X-3)/(8*(52*X-135))
  cosine=s.sqrt(X+s.Rational(1,4))*b/(2*s.sqrt((X-s.Rational(63,4))*a*c))
  q=simp(cosine+s.I*s.sqrt(1-cosine**2))
  parts=[{s.Rational(3,2):s.sqrt(s.Rational(5,8)),-s.Rational(5,2):-s.sqrt(s.Rational(3,8))},
         {s.Rational(5,2):-s.sqrt(s.Rational(3,8)),-s.Rational(3,2):-s.sqrt(s.Rational(5,8))},
         {s.Rational(9,2):s.sqrt(a)*q,s.Rational(1,2):s.sqrt(b),-s.Rational(7,2):s.sqrt(c)},
         {s.Rational(7,2):s.sqrt(c),-s.Rational(1,2):-s.sqrt(b),-s.Rational(9,2):s.sqrt(a)*q}]
  first_x=s.Rational(13,2);kind='central-four'
 basis=[vector(part)for part in parts]
 xs=[];x=first_x
 while x+1<=j:
  xs.append(x);va=s.sqrt((x+1)/(2*x+1));vb=s.sqrt(x/(2*x+1))
  basis.extend([vector({x:va,-x-1:vb}),vector({x+1:vb,-x:-va})]);x+=3
 assert len(basis)==n//3+1
 return {'id':f'{kind}-n{n}-k{len(basis)}-t1','n':n,'dimension':len(basis),'order':1,
  'label':'Central cyclic block with reflected outer pairs','basis':[coeff_record(v)for v in basis],
  'construction':{'type':kind,'centralProbabilities':{'a':str(a),'b':str(b),'c':str(c)},'phase':str(q),
    'magneticPairLowerValues':list(map(str,xs)),
    'argument':'The central block obeys exact zero compressed Jz and J± identities with a unit phase. Reflected outer pairs are separated by one unused magnetic level from the central block and from one another.'},
  'verification':validate(n,basis),'provenance':'General central-block construction; independently verified here by every reduced cross matrix.'}

def verify_general_blocks():
 X=s.Symbol('X',positive=True);Y=s.Symbol('Y',nonnegative=True);records=[]
 for kind in ['central-three','central-four']:
  if kind=='central-three':
   a=8*X/(9*(4*X-7));b=7*(4*X-15)/(18*(4*X-7));c=7*(4*X-3)/(18*(4*X-7));minimum=s.Rational(63,4)
   C2=(X+s.Rational(1,4))*b*b/(4*(X-s.Rational(35,4))*a*c)
   equations=[a+b+c-1,7*a+b-5*c,b*(X-s.Rational(3,4))-c*(X-s.Rational(15,4))]
  else:
   a=16*X/(52*X-135);b=27*(4*X-35)/(8*(52*X-135));c=45*(4*X-3)/(8*(52*X-135));minimum=s.Rational(99,4)
   C2=(X+s.Rational(1,4))*b*b/(4*(X-s.Rational(63,4))*a*c)
   equations=[a+b+c-1,9*a+b-7*c,5*b*(X-s.Rational(3,4))-3*c*(X-s.Rational(35,4))]
  assert all(s.factor(eq)==0 for eq in equations)
  sine2=s.factor(1-C2);num,den=s.fraction(sine2)
  shifted_num=s.Poly(s.expand(num.subs(X,Y+minimum)),Y);shifted_den=s.Poly(s.expand(den.subs(X,Y+minimum)),Y)
  assert all(c>0 for c in shifted_num.all_coeffs()) and all(c>0 for c in shifted_den.all_coeffs())
  for probability in [a,b,c]:
   numerator,denominator=s.fraction(probability)
   assert all(c>0 for c in s.Poly(s.expand(numerator.subs(X,Y+minimum)),Y).all_coeffs())
   assert all(c>0 for c in s.Poly(s.expand(denominator.subs(X,Y+minimum)),Y).all_coeffs())
  records.append({'type':kind,'minimumX':str(minimum),'probabilities':dict(zip(['a','b','c'],map(str,[a,b,c]))),
      'normalizationMeanAndCrossIdentitiesVerified':True,'phaseCosineSquared':str(s.factor(C2)),'phaseSineSquared':str(sine2),
      'shiftedPositiveNumerator':str(shifted_num.as_expr()),'shiftedPositiveDenominator':str(shifted_den.as_expr()),
      'allShiftedCoefficientsPositive':True,
      'finalPhaseEquation':'For central-three use negative cosine; for central-four use positive cosine. The remaining compressed J+ equation then vanishes by its defining phase formula.'})
 return records

def main():
 rows=[reflection_pairs(n)for n in list(range(2,21,2))+[5,11,17,22,23,24,26,35,42]]+[n7()]+[central_odd(n)for n in [9,13,15,19,21,25]]
 rows.sort(key=lambda r:(r['n'],r['id']))
 for row in rows:print(row['id'],'verified',flush=True)
 data={'schemaVersion':1,'subspaces':rows,
  'generalFamily':{'domain':'Every n>=4, as well as n=2.',
     'dimension':'floor(n/3)+1','order':1,
     'symbolicCentralBlockVerification':verify_general_blocks(),
     'remainingOddSizes':'For n=6q+1 (q>=1), use the central-three block and outer modules x=11/2+3r. For n=6q+3 (q>=1), use the central-four block and outer modules x=13/2+3r.',
     'integerSpinPrescription':'For j=3q or3q+1, include m=0 and modules x=2,5,...,3q-1. For j=3q+2, use modules x=1,4,...,3q+1 and omit m=0.',
     'halfIntegerSpinPrescription':'For n=6q+5, use modules x=3/2,9/2,...,3q+3/2 and omit m=0.',
     'scope':'This proves existence. Paper maximal-dimension labels remain conjectural unless a separate upper bound matches the dimension.'},
  'verification':'Every entry is checked symbolically for orthonormality and all k² reduced cross operators at t=1; a t=2 failure is recorded for each whole subspace.'}
 (ROOT/'t1-bases.json').write_text(json.dumps(data,indent=2)+'\n')
 print('Verified',len(rows),'exact subspaces.')
if __name__=='__main__':main()
