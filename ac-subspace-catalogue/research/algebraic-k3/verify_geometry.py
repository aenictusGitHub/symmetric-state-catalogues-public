#!/usr/bin/env python3
"""Certified Majorana multiplicities without algebraic primitive-element fields.
Run python3 verify_geometry.py --write to generate geometry and root disks.
Run python3 verify_geometry.py to independently check the stored exact rational
interval/Rouche certificates and the displayed Majorana coordinates.

Only generation needs mpmath. Proof arithmetic uses fractions.Fraction and
integer isqrt. A disk contains exactly one reduced-polynomial root (counting
multiplicity) by Rouche; disjoint disks exhaust the polynomial degree. Its
nonzero roots therefore are simple. Zero/infinity multiplicities follow from
exact coefficient support. Named algebraic parameters use the enclosures
proved by the corresponding subspace's separate exact existence certificate.
"""
from pathlib import Path
from fractions import Fraction as Q
from functools import reduce
from math import comb,gcd,isqrt,pi,sqrt
import cmath,hashlib,json,sys
import numpy as np
import sympy as s
BASE=Path(__file__).resolve().parent
SQRT_DIGITS=90
ZERO=(Q(0),Q(0));ONE=(Q(1),Q(1))
def radd(a,b):return a[0]+b[0],a[1]+b[1]
def rneg(a):return -a[1],-a[0]
def rmul(a,b):
 v=[x*y for x in a for y in b];return min(v),max(v)
def rsquare(a):return (Q(0) if a[0]<=0<=a[1] else min(a[0]**2,a[1]**2),max(a[0]**2,a[1]**2))
def rinv(a):
 assert not a[0]<=0<=a[1],('division by an interval containing zero',a)
 return 1/a[1],1/a[0]
def rsqrt(a):
 assert a[0]>=0
 scale=10**SQRT_DIGITS
 def floorroot(q):return isqrt(q.numerator*scale*scale//q.denominator)
 lo=floorroot(a[0]);hi=floorroot(a[1]);upper=Q(hi,scale)
 if upper*upper<a[1]:upper=Q(hi+1,scale)
 return Q(lo,scale),upper
def cpoint(re,im=0):return ((Q(re),Q(re)),(Q(im),Q(im)))
CZERO=cpoint(0);CONE=cpoint(1)
def cadd(a,b):return radd(a[0],b[0]),radd(a[1],b[1])
def cmul(a,b):return radd(rmul(a[0],b[0]),rneg(rmul(a[1],b[1]))),radd(rmul(a[0],b[1]),rmul(a[1],b[0]))
def cinv(a):
 den=rinv(radd(rsquare(a[0]),rsquare(a[1])));return rmul(a[0],den),rmul(rneg(a[1]),den)
def cpow(a,k):
 if k<0:return cpow(cinv(a),-k)
 out=CONE
 while k:
  if k&1:out=cmul(out,a)
  a=cmul(a,a);k//=2
 return out
def absupper(a):return max(abs(v) for v in a[0])+max(abs(v) for v in a[1])
def abslower(a):
 def component(b):return 0 if b[0]<=0<=b[1] else min(abs(b[0]),abs(b[1]))
 return max(component(a[0]),component(a[1]))
def interval_eval(expression,parameters):
 cache={}
 def go(e):
  if e in cache:return cache[e]
  if e==s.I:out=cpoint(0,1)
  elif e.is_Rational:out=cpoint(Q(int(e.p),int(e.q)))
  elif e.is_Symbol:
   assert str(e) in parameters,('missing parameter',e)
   bounds=tuple(map(Q,parameters[str(e)]));assert bounds[0]<=bounds[1];out=(bounds,ZERO)
  elif e.is_Add:out=reduce(cadd,(go(a) for a in e.args),CZERO)
  elif e.is_Mul:out=reduce(cmul,(go(a) for a in e.args),CONE)
  elif e.is_Pow:
   b,p=e.args;v=go(b)
   if p.is_Integer:out=cpow(v,int(p))
   elif p.is_Rational and p.q==2:
    assert v[1]==ZERO
    out=cpow((rsqrt(v[0]),ZERO),int(p.p))
   else:raise ValueError(('unsupported exponent',e))
  else:raise ValueError(('unsupported expression',e))
  cache[e]=out;return out
 return go(expression)
def source(n):
 data=json.loads((BASE/f'n{n}'/'basis.json').read_text());return data.get('subspaces',[data])[0]
def basis_polynomial(space,row):
 n=space['n'];parameters=space.get('algebraicDefinition',{}).get('parameterIntervals',{})
 locals_={k:s.Symbol(k,real=True) for k in parameters}
 expressions=[s.sympify(e,locals=locals_,rational=True) for e in row['exactCoefficients']]
 coeffs=[interval_eval((-1)**r*s.sqrt(s.binomial(n,r))*e,parameters) for r,e in enumerate(expressions)]
 support=[r for r,v in enumerate(coeffs) if v!=CZERO];lo,hi=min(support),max(support)
 assert abslower(coeffs[lo])>0 and abslower(coeffs[hi])>0
 step=reduce(gcd,(r-lo for r in support),0) or 1
 assert all(coeffs[r]==CZERO for r in range(n+1) if (r-lo)%step)
 reduced=[coeffs[r] for r in range(lo,hi+1,step)]
 return expressions,parameters,reduced,lo,hi,step

def taylor(poly,center):
 degree=len(poly)-1
 return [reduce(cadd,(cmul(cpoint(comb(k,j)),cmul(poly[k],cpow(center,k-j))) for k in range(j,degree+1)),CZERO) for j in range(degree+1)]
def rounded_up(q,digits=50):
 scale=10**digits;return Q((q.numerator*scale+q.denominator-1)//q.denominator,scale)
def rounded_down(q,digits=50):
 scale=10**digits;return Q(q.numerator*scale//q.denominator,scale)
def disk_check(poly,disk):
 center=cpoint(Q(disk['center'][0]),Q(disk['center'][1]));radius=Q(disk['radius']);assert radius>0
 coeff=taylor(poly,center)
 left=absupper(coeff[0])+sum(absupper(coeff[j])*radius**j for j in range(2,len(coeff)))
 right=abslower(coeff[1])*radius
 assert left<right,('Rouche inequality failed',float(left),float(right))
 assert abslower(center)>radius,'disk includes zero'
 return left,right

def generate_disks(poly):
 import mpmath as mp
 mp.mp.dps=110
 def mpq(v):return mp.mpf(v.numerator)/v.denominator
 center=[mp.mpc(mpq(sum(c[0],Q(0))/2),mpq(sum(c[1],Q(0))/2)) for c in poly]
 roots=mp.polyroots(list(reversed(center)),maxsteps=3000,error=False,extraprec=160)
 roots.sort(key=lambda z:(float(z.real),float(z.imag)));out=[]
 scale=10**70
 for z in roots:
  disk={'center':[str(Q(int(mp.nint(z.real*scale)),scale)),str(Q(int(mp.nint(z.imag*scale)),scale))],'radius':str(Q(1,10**15))}
  left,right=disk_check(poly,disk)
  disk.update(roucheLeftUpper=str(rounded_up(left)),roucheRightLower=str(rounded_down(right)))
  out.append(disk)
 return out

def point(root,m):
 if root is None:return {'x':0.,'y':0.,'z':1.,'multiplicity':m}
 q=1+abs(root)**2;return {'x':2*root.real/q,'y':-2*root.imag/q,'z':(abs(root)**2-1)/q,'multiplicity':m}
def points_from_disks(disks,n,lo,hi,step):
 points=[]
 if lo:points.append(point(0j,lo))
 if hi<n:points.append(point(None,n-hi))
 for disk in disks:
  root=complex(float(Q(disk['center'][0])),float(Q(disk['center'][1])))
  for k in range(step):points.append(point(cmath.rect(abs(root)**(1/step),(cmath.phase(root)+2*pi*k)/step),1))
 assert sum(p['multiplicity'] for p in points)==n
 return points

def reconstruct(points,n):
 poly=np.array([1.+0j])
 for p in points:
  a=sqrt(max(0.,(1+p['z'])/2));b=(p['x']+1j*p['y'])/(2*a) if a>1e-10 else 1.+0j
  for _ in range(p['multiplicity']):poly=np.convolve(poly,[a,-b])
 return np.array([(-1)**r*poly[r]/sqrt(comb(n,r)) for r in range(n+1)])
def ray_error(a,b):
 a=np.asarray(a,complex);b=np.asarray(b,complex);a/=np.linalg.norm(a);b/=np.linalg.norm(b);overlap=np.vdot(a,b)
 return float(np.linalg.norm(a-b*np.conj(overlap)/abs(overlap)))
def hashes(row,parameters):
 return hashlib.sha256(json.dumps({'exactCoefficients':row['exactCoefficients'],'parameterIntervals':parameters},sort_keys=True,separators=(',',':')).encode()).hexdigest()

def compare_rendering(actual,expected,path=()):
 """Exact certificate comparison; only explicitly listed rendering floats may drift.
 NumPy/SymPy evaluation in a shared process can change the final rounding bit.
 This tolerance never applies to intervals, root disks, hashes or multiplicities.
 """
 if isinstance(actual,dict):
  assert isinstance(expected,dict) and actual.keys()==expected.keys(),path
  for key in actual:compare_rendering(actual[key],expected[key],path+(key,))
 elif isinstance(actual,list):
  assert isinstance(expected,list) and len(actual)==len(expected),path
  for i,(a,b) in enumerate(zip(actual,expected)):compare_rendering(a,b,path+(i,))
 else:
  point_coordinate=len(path)>=3 and path[-3]=='points' and path[-1] in ['x','y','z']
  diagnostic=len(path)>=2 and path[-2]=='diagnostics' and path[-1] in ['projectiveDickeReconstructionError','exactFormulaVsSavedCoefficientsError','maxUnitSphereError','minimumDistinctPointDistance']
  summary=len(path)>=2 and path[-2]=='verification' and path[-1] in ['maximumProjectiveDickeReconstructionError','minimumDistinctPointDistance']
  if (point_coordinate or diagnostic or summary) and isinstance(actual,float) and isinstance(expected,float):
   assert abs(actual-expected)<=1e-14,('rendering differs beyond rounding tolerance',path,actual,expected)
  else:assert actual==expected,('exact certificate data differs',path,actual,expected)

def run(write=False,cases=(14,15,17),source_fn=source,output_dir=BASE):
 saved=None if write else json.loads((output_dir/'geometry.json').read_text());records={};maxerr=0.;mindistance=2.;diskcount=0
 for n in cases:
  space=source_fn(n)
  for index,row in enumerate(space['basis']):
   id_=space['id']+'-basis-'+str(index+1)
   expressions,parameters,poly,lo,hi,step=basis_polynomial(space,row)
   disks=generate_disks(poly) if write else saved['records'][id_]['rootDisks']
   assert len(disks)==len(poly)-1
   for i,disk in enumerate(disks):
    left,right=disk_check(poly,disk)
    assert left<=Q(disk['roucheLeftUpper'])<Q(disk['roucheRightLower'])<=right
    for other in disks[:i]:
     a=[Q(v) for v in disk['center']];b=[Q(v) for v in other['center']]
     assert sum((a[k]-b[k])**2 for k in range(2))>(Q(disk['radius'])+Q(other['radius']))**2
   points=points_from_disks(disks,n,lo,hi,step);rendered=[complex(*c) for c in row['coefficients']]
   restored=reconstruct(points,n);error=ray_error(rendered,restored);assert error<1e-10,(id_,error)
   maxerr=max(maxerr,error);minimum=min(sqrt(sum((a[k]-b[k])**2 for k in ['x','y','z'])) for i,a in enumerate(points) for b in points[:i]);mindistance=min(mindistance,minimum)
   assert minimum>1e-7;uniterror=max(abs(sum(p[k]**2 for k in ['x','y','z'])-1) for p in points)
   # Ensure saved rendering agrees with the exact expressions' enclosure center.
   centers={s.Symbol(k,real=True):s.Rational(Q(v[0])+Q(v[1]))/2 for k,v in parameters.items()}
   exactrender=[complex(s.N(e.subs(centers),35)) for e in expressions];formulaerror=ray_error(rendered,exactrender);assert formulaerror<1e-12
   record={'id':id_,'n':n,'label':row.get('label',f'Basis vector {index+1}'),'points':points,'partition':sorted([p['multiplicity'] for p in points],reverse=True),'sourceExactCoefficients':row['exactCoefficients'],'sourceParameterIntervals':parameters,'sourceSHA256':hashes(row,parameters),'rootDisks':disks,'coordinateMethod':'Exact rational interval arithmetic and Rouche disks certify all nonzero roots simple; numerical disk centers render the sphere.','multiplicityEvidence':'Exact rational Rouche certificate and exact pole support.','evidence':'Exact AC subspace basis vector; coordinates are a numerical rendering of its algebraic coefficients.','diagnostics':{'finiteRootsCertifiedSimple':True,'multiplicityCertifiedExact':True,'squareFreeFactorizationExact':False,'rootCertificationMethod':'rational interval Rouche theorem','zeroMultiplicity':lo,'infinityMultiplicity':n-hi,'nonzeroRootStep':step,'reducedPolynomialDegree':len(poly)-1,'disjointCertifiedRootDisks':len(disks),'roucheDiskRadius':disks[0]['radius'] if disks else None,'intervalSquareRootDigits':SQRT_DIGITS,'rootPrecisionDecimalDigits':70,'projectiveDickeReconstructionError':error,'exactFormulaVsSavedCoefficientsError':formulaerror,'maxUnitSphereError':uniterror,'minimumDistinctPointDistance':minimum}}
   if not write:
    compare_rendering(record,saved['records'][id_],('records',id_))
   records[id_]=record;diskcount+=len(disks);print(id_,error,flush=True)
 report={'passed':True,'recordCount':len(records),'certifiedReducedRootDisks':diskcount,'allFiniteNonzeroRootsCertifiedSimple':True,'allMultiplicitiesCertified':True,'maximumProjectiveDickeReconstructionError':maxerr,'minimumDistinctPointDistance':mindistance,'certificateArithmetic':'Exact rational fractions; square roots enclosed using integer isqrt; no floating-point inequalities used for multiplicity proofs.'}
 result={'schemaVersion':1,'convention':{'polynomial':'P(z)=sum_r (-1)^r sqrt(binomial(n,r))*c_r*z^r','sphere':'(2Re(z), -2Im(z), |z|^2-1)/(1+|z|^2); infinity is north.'},'records':records,'verification':report}
 if write:(output_dir/'geometry.json').write_text(json.dumps(result,indent=2)+'\n')
 else:compare_rendering(result,saved)
 (output_dir/'geometry-verification.json').write_text(json.dumps(report,indent=2)+'\n')
 return report

def verify():return run(False)
if __name__=='__main__':print(json.dumps(run('--write' in sys.argv),indent=2))
