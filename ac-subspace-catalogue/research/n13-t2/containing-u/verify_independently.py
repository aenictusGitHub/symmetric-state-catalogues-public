#!/usr/bin/env python3
"""Independent exact physical-matrix verification. Requires SymPy."""
import sympy as s,json,itertools,math,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
start=time.time();d=json.loads((ROOT/'certificate.json').read_text());n=13;j=s.Rational(13,2)
u=s.zeros(14,1);u[0]=s.sqrt(s.Rational(5,18));u[9]=s.sqrt(s.Rational(13,18))
Z=s.diag(*[j-r for r in range(14)]);L=s.zeros(14)
for r in range(13):L[r+1,r]=s.sqrt((r+1)*(13-r))
R=L.T;ops=[Z,Z*Z-s.Rational(65,4)*s.eye(14),L,R,Z*L+L*Z,Z*R+R*Z,L*L,R*R]
images=s.Matrix.hstack(u,*[A*u for A in ops]);assert images.rank()==7
H=s.Matrix.hstack(*[s.Matrix([s.sympify(v)for v in column])for column in d['kernelBasisDicke']]);G=(H.H*H).applyfunc(s.simplify)
assert H.rank()==7 and G==s.diag(*d['kernelGram'])
assert (images.H*H).applyfunc(s.simplify)==s.zeros(9,7)
# Explicit charge audit: supports of each h_a have one residue modulo 9.
charges=[]
for a in range(7):
 support=[r for r in range(14)if H[r,a]!=0];q=support[0]%9;assert all(r%9==q for r in support);charges.append(q)
assert charges==d['charges']
A=[(H.H*T*H).applyfunc(s.simplify)for T in ops];assert A==[s.Matrix([[s.sympify(v)for v in row]for row in T])for T in d['compressedOperators']]
pairs=list(itertools.combinations(range(7),2));assert [list(p)for p in pairs]==d['pairs'];sectors={}
for a,(i,j)in enumerate(pairs):sectors.setdefault((charges[i]+charges[j])%9,[]).append(a)
sectors=dict(sorted(sectors.items()));assert sectors=={int(k):v for k,v in d['sectors'].items()}
variables=[(q,a,b)for q,ind in sectors.items()for a in range(len(ind))for b in range(a,len(ind))];assert [list(v)for v in variables]==d['variables']
W=[]
for i,j in pairs:
 M=s.zeros(7);M[i,j]=1;M[j,i]=-1;W.append(M)
columns=[]
for q,a,b in variables:
 aa,bb=sectors[q][a],sectors[q][b];column=[]
 for T in A:
  M=W[aa].T*T.T*W[bb]
  if a!=b:M+=W[bb].T*T.T*W[aa]
  column.extend(M)
 columns.append(column)
M=s.Matrix.hstack(*[s.Matrix(c)for c in columns]);rows=[];metadata=[]
for pos in range(M.rows):
 row=list(M.row(pos))
 if not any(row):continue
 den=math.lcm(*[int(x.q)for x in row]);row=[int(x*den)for x in row];g=math.gcd(*row);rows.append([v//g for v in row]);metadata.append([pos//49,pos%49//7,pos%7])
assert rows==d['allPhysicalRows']and metadata==d['allRowMetadata']
selected=[rows[i]for i in d['selectedIndices']];assert selected==d['constraintRows']
B=s.Matrix(selected);K=s.Matrix(d['kernel']);assert B.rank()==40 and K.rank()==8 and B*K==s.zeros(40,8)
assert s.Matrix(rows)*K==s.zeros(len(rows),8)
diag=[i for i,(_,a,b)in enumerate(variables)if a==b];assert diag==d['diagonalVariableIndices'];c=d['positiveCoefficients'];assert len(c)==21 and min(c)>0 and s.Matrix([c])*K[diag,:]==s.zeros(1,8)
report={'passed':True,'crossImageRank':7,'fixedStateOrthogonalComplementDimension':7,'compressedOperatorsExactlyVerified':8,'C9ChargesExactlyVerified':charges,'physicalRowsIndependentlyReconstructed':len(rows),'invariantCoordinates':len(variables),'rank':40,'kernelDimension':8,'strictPositiveWeights':21,'scope':'No 3D second-order AC subspace containing the specified u, or a collective rotation thereof; not a global N13,t2 exclusion.','runtimeSeconds':time.time()-start}
(ROOT/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
