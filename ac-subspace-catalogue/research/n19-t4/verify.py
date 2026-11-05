#!/usr/bin/env python3
"""Nineteen-qubit AC4 witness: independent exact standard-library arithmetic."""
from collections import defaultdict
from fractions import Fraction as F
from math import comb,gcd,isqrt,prod
from pathlib import Path
import ast,hashlib,json
HERE=Path(__file__).resolve().parent
# A term (epsilon,radicand) represents i^epsilon sqrt(radicand).
# All coefficient radicands below are square free; their factorization is checked.
FACTORS={194:[2,97],17290:[2,5,7,13,19],2918881141085:[5,7,13,19,109,431,7187],266:[2,7,19],7410:[2,3,5,13,19],323323:[7,11,13,17,19],218332309353158:[2,7,11,13,17,19,109,431,7187]}
C={0:{(0,194):F(1,40)},6:{(0,17290):F(-9673,5605600),(1,2918881141085):F(-1,2802800)},10:{(0,266):F(1,160)},14:{(0,7410):F(1,160)},18:{(0,323323):F(-2837,10210200),(1,218332309353158):F(-1,40840800)}}
def add(a,b):
 c=defaultdict(F,a)
 for key,v in b.items():c[key]+=v
 return {key:v for key,v in c.items()if v}
def mul(a,b):
 c=defaultdict(F)
 for (e,r),v in a.items():
  for (f,s),w in b.items():
   g=gcd(r,s);c[((e+f)%2,(r//g)*(s//g))]+=v*w*g*(-1 if e+f==2 else 1)
 return {key:v for key,v in c.items()if v}
def conj(a):return {key:v*(-1 if key[0] else 1)for key,v in a.items()}
def sqrt_rational(v):
 v=F(v);assert v>=0
 if not v:return {}
 a=v.numerator*v.denominator;factor=1;rad=1;p=2
 while p*p<=a:
  power=0
  while a%p==0:a//=p;power+=1
  factor*=p**(power//2)
  if power%2:rad*=p
  p+=1
 if a>1:rad*=a
 return {(0,rad):F(factor,v.denominator)}
def parse_coefficient(text):
 def run(node):
  if isinstance(node,ast.Constant)and isinstance(node.value,int):return {(0,1):F(node.value)}if node.value else {}
  if isinstance(node,ast.Name)and node.id=='I':return {(1,1):F(1)}
  if isinstance(node,ast.UnaryOp)and isinstance(node.op,ast.USub):return {k:-v for k,v in run(node.operand).items()}
  if isinstance(node,ast.Call)and isinstance(node.func,ast.Name)and node.func.id=='sqrt'and len(node.args)==1 and not node.keywords:
   value=run(node.args[0]);assert list(value)==[(0,1)]
   q=value[0,1];assert q.denominator==1 and q.numerator in FACTORS
   return {(0,q.numerator):F(1)}
  if isinstance(node,ast.BinOp):
   a,b=run(node.left),run(node.right)
   if isinstance(node.op,ast.Add):return add(a,b)
   if isinstance(node.op,ast.Sub):return add(a,{k:-v for k,v in b.items()})
   if isinstance(node.op,ast.Mult):return mul(a,b)
   if isinstance(node.op,ast.Div):
    assert list(b)==[(0,1)]and b[0,1]
    return {k:v/b[0,1]for k,v in a.items()}
  raise AssertionError('Unrecognized coefficient expression')
 return run(ast.parse(text,mode='eval').body)
def rho(t):
 out=[]
 for a in range(t+1):
  row=[]
  for b in range(t+1):
   value={}
   for l in range(20-t):
    r,s=a+l,b+l
    if r not in C or s not in C:continue
    scale=F(comb(t,a)*comb(t,b)*comb(19-t,l)**2,comb(19,r)*comb(19,s))
    value=add(value,mul(mul(C[r],conj(C[s])),sqrt_rational(scale)))
   row.append(value)
  out.append(row)
 return out
def verify(directory=None):
 root=Path(directory or HERE)
 for rad,ps in FACTORS.items():
  assert len(set(ps))==len(ps)and prod(ps)==rad
  assert all(p>=2 and all(p%d for d in range(2,isqrt(p)+1))for p in ps)
 expected={0:F(97,800),6:F(1083,2560),10:F(133,12800),14:F(741,2560),18:F(399,2560)}
 for r,c in C.items():assert mul(c,conj(c))=={(0,1):expected[r]}
 assert sum(expected.values())==1
 moments=[]
 for power in range(5):
  value=sum(w*r**power for r,w in expected.items());target=F(sum(i**power for i in range(20)),20);assert value==target
  moments.append({'power':power,'value':str(value)})
 checks=[]
 for t in range(1,5):
  matrix=rho(t)
  for a in range(t+1):
   for b in range(t+1):assert matrix[a][b]==({(0,1):F(1,t+1)}if a==b else {}),(t,a,b,matrix[a][b])
  checks.append({'t':t,'identityVerified':True,'entriesChecked':(t+1)**2,'frobeniusDistanceSquared':'0'})
 fifth=rho(5);distance={}
 for a in range(6):
  for b in range(6):
   delta=add(fifth[a][b],{(0,1):F(-1,6)}if a==b else {})
   distance=add(distance,mul(delta,conj(delta)))
 assert list(distance)==[(0,1)]and distance[0,1]>0
 # t=5 is already obstructed by the diagonal moment for this particular state.
 moment5=sum(w*r**5 for r,w in expected.items())-F(sum(i**5 for i in range(20)),20);assert moment5
 result={'passed':True,'n':19,'t':4,'exactAchievedOrder':4,'probabilities':{str(k):str(v)for k,v in expected.items()},'exactDiagonalMoments':moments,'reductions':checks,'exactReducedMatrixEntriesChecked':sum(c['entriesChecked']for c in checks),'fifthOrderFrobeniusDistanceSquared':str(distance[0,1]),'fifthRawMomentDifference':str(moment5),'method':'Independent full normalized-Dicke splitting; exact fractions and square-free radicals, standard library only. All reduced matrix entries through order5 computed.','scope':'Exact existence of a fourth-order AC line. Does not assert that AC4 is the highest order possible on19 qubits, or that its one-dimensional span is the largest AC4 subspace.'}
 if(root/'basis.json').exists():
  raw=(root/'basis.json').read_bytes();data=json.loads(raw);coeffs=data['basis'][0]['exactCoefficients'];assert len(coeffs)==20
  assert all(parse_coefficient(c)==C.get(r,{})for r,c in enumerate(coeffs))
  result['exactSourceCoefficientIdentitiesChecked']=20;result['sourceBasisSHA256']=hashlib.sha256(raw).hexdigest()
 return result
if __name__=='__main__':
 report=verify();(HERE/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
