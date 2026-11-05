"""Complete the rigorous identity and positivity verification."""
from exact_field import *
from flint import arb,arb_mat,ctx,fmpq_mat,fmpq_poly
from fractions import Fraction as F
from math import comb
from collections import defaultdict
import json,hashlib,time
ROOT=Path(__file__).parent;ctx.prec=512;start=time.time();cert=json.load(open(ROOT/'sos-certificate.json'));check=json.load(open(ROOT/'exact-checks.json'));assert all(check[k]for k in ['covarianceVerified','kernelVanishingVerified','rightInvarianceVerified','exactAttainmentVerified','exactStationarityVerified']);assert check['invariantDimension']==290;assert cert['gramBase']=='identity'and cert['pivotInterpretation']=='additive symmetric corrections'
lo,hi=map(fmpq,DATA['rootInterval']);assert MOD(lo)*MOD(hi)<0
for _ in range(1150):
 mid=(lo+hi)/2
 if MOD(lo)*MOD(mid)<0:hi=mid
 else:lo=mid
ball=arb(lo).union(arb(hi))
def ev(z):
 v=arb(0)
 for c in reversed(K(z).p.coeffs()):v=v*ball+arb(c)
 return v
H=json.load(open(ROOT/'basis-evaluations.json'));J=json.load(open(ROOT/'exact-jets.json'));assert len(H)==286 and all(len(row)==290 for row in H)and len(J)==290
C=arb_mat([[arb(fmpq(v))for v in row]for row in H]+[[ev(J[i][j])for i in range(290)]for j in range(4)]);det=C.det();print('INJECTIVITY DETERMINANT',det,flush=True);assert not det.contains(0)
values=json.load(open(ROOT/'interval-evaluations.json'));assert len(values)==286 and all(len(row)==287 for row in values);A=arb_mat([[arb(v)for v in row[1:]]for row in values]);poly=[(tuple(m),int(c))for m,c in json.load(open(ROOT/'target-polynomial.json'))]
def mul(a,b):return(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
target=[]
for i,point in enumerate(cert['evaluationPoints']):
 xx=[tuple(x)for x in point]+[(a,-b)for a,b in point];q=sum(comb(11,j)*(a*a+b*b)for j,(a,b)in enumerate(point));total=[0,0]
 for mon,c in poly:
  v=(1,0)
  for j in mon:v=mul(v,xx[j])
  total[0]+=c*v[0];total[1]+=c*v[1]
 assert total[1]==0;target.append(arb(fmpq(total[0],q*q))-ev(BOUND)-arb(values[i][0]))
solution=A.solve(arb_mat([[x]for x in target]),algorithm='precond');print('PASS exact linear system has a unique rigorously enclosed solution',time.time()-start,flush=True)
dfn=json.load(open(ROOT/'covariants.json'));variables=[(ib,a,b)for ib,B0 in enumerate(dfn['blocks'])for b in range(B0['dimension'])for a in range(b+1)];piv=cert['pivotVariables'];assert len(piv)==len(set(piv))==286 and all(0<=i<len(variables)for i in piv);G=[[[arb(a==b)for b in range(B0['dimension'])]for a in range(B0['dimension'])]for B0 in dfn['blocks']]
for j,v in enumerate(piv):
 ib,a,b=variables[v];G[ib][a][b]+=solution[j,0]
 if a!=b:G[ib][b][a]+=solution[j,0]
bounds=[]
for ib,H0 in enumerate(G):
 for a in range(len(H0)):
  lower=H0[a][a]-sum((abs(H0[a][b])for b in range(len(H0))if b!=a),arb(0));assert lower>0;bounds.append(lower)
 print('PASS strictly diagonally dominant positive Gram block',ib,flush=True)
assert all(b>arb(fmpq(99,100))for b in bounds);minimum='99/100';print('PROVED GLOBAL MINIMUM',ev(BOUND));print('minimum positivity lower bound',minimum,flush=True)
value=json.load(open(ROOT/'value-definition.json'));vp=fmpq_poly([fmpq(v)for v in value['minimalPolynomialCoefficientsDescending'][::-1]]);assert vp.degree()==11
powers=[K(1)]
for j in range(11):powers.append(powers[-1]*BOUND)
assert not sum((K(vp[j])*powers[j]for j in range(12)),K());assert fmpq_mat([[powers[j].p[i]for j in range(11)]for i in range(11)]).det()!=0
vlow,vhigh=map(fmpq,value['rootInterval']);assert ev(BOUND)>arb(vlow)and ev(BOUND)<arb(vhigh);assert vp(vlow)*vp(vhigh)<0
vi=arb(vlow).union(arb(vhigh));dv=arb(0)
for c in reversed(vp.derivative().coeffs()):dv=dv*vi+arb(c)
assert not dv.contains(0)
print('PASS exact degree-eleven minimum polynomial and isolating interval',flush=True)
report=dict(n=11,sector=6,direction='min',degree=10,globalOptimalityProven=True,fullPolynomialIdentityVerified=True,attainmentVerifiedExact=True,exactStationarityVerified=True,rightSymmetry='Kramers SU(2) mixing; exact right-generator invariance checked',invariantDimension=290,coefficientConstraints=286,vanishingConstraints=4,positiveDefiniteGramBlocks=len(G),certifiedGershgorinLowerBound=minimum,arithmetic='Exact rational/algebraic checks; rigorously bounded integer evaluation; directed-rounding FLINT/Arb arithmetic',floatingPointToleranceUsed=False,fieldEvaluationPrecisionBits=1024,polynomialEvaluationPrecisionBits=256,linearSolvePrecisionBits=512,boundEnclosure=ev(BOUND).str(100),injectivityDeterminant=str(det),identity='The 286 evaluation equations and four exact stationary vanishing conditions are injective on the complete 290-dimensional jointly invariant space.',gramDefinition='Each Gram matrix is the identity plus the symmetric pivot corrections uniquely defined by the exact evaluation equations.',inputSha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()for name in ['source-basis.json','source-kernels.json','holomorphic-basis.json','covariants.json','sos-certificate.json','state-parameters.json','exact_field.py','exact_checks.py','integer_evaluations.py','basis_values.py','interval_evaluations.py','arb_helper.c','verify_final.py','check_couplings.py','verify_certificate.py','requirements.txt','value-definition.json']})
(ROOT/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
(ROOT/'pivot-enclosures.json').write_text(json.dumps([solution[i,0].str(90)for i in range(286)],indent=2)+'\n')
