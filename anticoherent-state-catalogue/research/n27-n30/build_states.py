"""Prepare exact attainment records for N27--30.

Requires SymPy, mpmath and NumPy only for coefficient/coordinate renderings.
The saved proofs are checked by the independent scripts in each subfolder.
When installed as research/n27-n30, the default source is data/catalogue.json.
"""
from pathlib import Path
from fractions import Fraction as Q
from functools import reduce
from math import gcd,comb
import argparse,copy,json
import mpmath as mp
import numpy as np
import sympy as s
HERE=Path(__file__).resolve().parent

def geometry(record):
 n=record['n'];c=list(map(s.sympify,record['exactCoefficients']));support=[r for r in range(n+1) if c[r]!=0]
 lo,hi=min(support),max(support);step=reduce(gcd,[r-lo for r in support]) or 1;degree=(hi-lo)//step
 with mp.workdps(100):
  coeff=[mp.mpc(0)]*(degree+1)
  for r in support:
   value=s.N((-1)**r*s.sqrt(comb(n,r))*c[r],105)
   coeff[(r-lo)//step]=mp.mpc(str(s.re(value)),str(s.im(value)))
  coeff=[x/coeff[-1] for x in coeff]
  try:roots=mp.polyroots(list(reversed(coeff)),maxsteps=2000,extraprec=150)
  except (mp.libmp.libhyper.NoConvergence,ValueError,ZeroDivisionError):
   roots=[];deriv=[j*coeff[j] for j in range(1,len(coeff))]
   for z in np.roots(np.array([complex(x) for x in reversed(coeff)])):
    z=mp.mpc(complex(z))
    for _ in range(100):
     delta=mp.polyval(list(reversed(coeff)),z)/mp.polyval(list(reversed(deriv)),z);z-=delta
     if abs(delta)<mp.mpf('1e-90'):break
    roots.append(z)
  rebuilt=[mp.mpc(1)]
  for z in roots:
   nxt=[mp.mpc(0)]*(len(rebuilt)+1)
   for j,v in enumerate(rebuilt):nxt[j]-=z*v;nxt[j+1]+=v
   rebuilt=nxt
  error=max(abs(a-b) for a,b in zip(rebuilt,coeff))/max(abs(x) for x in coeff)
  assert error<mp.mpf('1e-70')
  points=[]
  if lo:points.append({'x':0.,'y':0.,'z':-1.,'multiplicity':lo})
  if hi<n:points.append({'x':0.,'y':0.,'z':1.,'multiplicity':n-hi})
  for z in roots:
   for k in range(step):
    root=mp.root(abs(z),step)*mp.exp(1j*(mp.arg(z)+2*mp.pi*k)/step);r2=abs(root)**2;den=1+r2
    points.append({'x':float(2*mp.re(root)/den),'y':float(-2*mp.im(root)/den),'z':float((r2-1)/den),'multiplicity':1})
  return points,str(error)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,default=HERE.parents[1]/'data/catalogue.json');args=parser.parse_args()
 states=copy.deepcopy(json.loads((HERE/'radical-witnesses/states.json').read_text())['states'])
 checks=json.loads((HERE/'radical-witnesses/independent-verification.json').read_text())['states']
 for rec in states:
  n,t=rec['n'],rec['achievedOrder'];proof=next(x for x in checks if x['n']==n)
  rec['stateLatex']=r'\begin{gathered}'+rec['stateLatex']+r'\\[6pt]'+rec['parametersLatex']+r'\end{gathered}'
  rec['majoranaPoints'],error=geometry(rec);rec['partition']=proof['majoranaPartition']
  assert sorted([p['multiplicity'] for p in rec['majoranaPoints']],reverse=True)==rec['partition']
  rec['majoranaCoordinatePolicy']='numeric-roots-retain-explicit-poles'
  rec['verification']['majoranaMultiplicityProof']='Exact reciprocal-polynomial / polynomial-gcd checks in the independent verifier.'
  rec['verification']['majoranaRootReconstructionRelativeError']=error
  rec['references']=[{'title':f'Exact analytical order-{t} state at N={n}','url':f'./research/n27-n30/index.html#n{n}'},{'title':'Radical constructions and exact reduction proof','url':'./research/n27-n30/radical-witnesses/proof.tex'},{'title':'Independent exact verification','url':'./research/n27-n30/radical-witnesses/independent-verification.json'}]
 original=next(x for x in json.loads(args.source.read_text())['states'] if x['id']=='liege-acs.28.6')
 cert=json.loads((HERE/'n28-t6/certificate.json').read_text());report=json.loads((HERE/'n28-t6/independent-verification.json').read_text());center=list(map(Q,cert['center']));mp.mp.dps=150
 approx=[]
 for r in range(29):
  factor=mp.sqrt(comb(28,r))/cert['scaleFactors'][r]
  approx.append([mp.nstr(factor*mp.mpf(center[k].numerator)/center[k].denominator,100) for k in (r,29+r)])
 n28={'id':'certified-algebraic-n28-t6','n':28,'achievedOrder':6,'witnessesOrders':list(range(1,7)),'label':'Certified algebraic state · order 6','status':'proven-existence','coefficientEvidence':'implicit-algebraic','stateLatex':r'|\psi\rangle=\sum_{r=0}^{28}\frac{\sqrt{\binom{28}{r}}}{d_r}\bigl(x_r^*+ix_{29+r}^*\bigr)|D_{28}^{(r)}\rangle','support':list(range(29)),'coefficients':[[float(x),float(y)] for x,y in approx],'highPrecisionCoefficients':approx,'majoranaPoints':copy.deepcopy(original['majoranaPoints']),'partition':copy.deepcopy(original['partition']),'majoranaCoordinatePolicy':'Numeric rendering copied from the nearby refined Liège state; the exact state is defined by the certified polynomial root.','normalization':'Exact rho_6 = I_7/7 implies sum_r |c_r|^2 = 1.','coefficientConvention':'Normalized Dicke basis; r counts excitations. Floating components render the rational certificate center.','implicitAlgebraicDefinition':{'certificateUrl':'./research/n27-n30/n28-t6/certificate.json','proofUrl':'./research/n27-n30/n28-t6/proof.tex','equations':'49 real rational quadratics equivalent to rho_6=I_7/7 after binomial scaling.','variableConvention':cert['variableConvention'],'activeVariables':cert['activeVariables'],'fixedVariables':cert['fixedVariables'],'radius':cert['radius'],'rootDefinition':'Unique real zero in the active-coordinate infinity-norm box around the saved rational center, with nine coordinates fixed to their rational center values.','centerCoefficientErrorBound':'2e-65 per real or imaginary Dicke component, relative to the exact rational-center conversion; displayed floating values have additional rounding.','closedFormRadicalsClaimed':False},'derivedFromStateId':'liege-acs.28.6','sourceLibraryRecordId':'acs.28.6','provenance':'Certified exact algebraic zero near refined Liège acs.28.6, with exact rational Banach contraction proof.','references':[{'title':'Certified sixth-order existence at N=28','url':'./research/n27-n30/index.html#n28'},{'title':'Exact algebraic existence proof','url':'./research/n27-n30/n28-t6/proof.tex'},{'title':'Rational root certificate','url':'./research/n27-n30/n28-t6/certificate.json'},{'title':'Independent rational verification','url':'./research/n27-n30/n28-t6/independent-verification.json'}],'notes':'Exactly sixth-order AC and not seventh-order AC, proved by a rational isolating-root certificate and a next-order interval. This is an exact implicit algebraic state; elementary radical components are not claimed. Global maximum order requires a separate universal exclusion.','verification':{'method':'Exact rational reduction equations, reconstructed Jacobian and strict Banach contraction bounds.','exactExistenceMethod':'rational-contraction','exactNormalization':True,'exactOrderVerified':True,'globalMaximumOrderProven':False,'existenceCertificate':'./research/n27-n30/n28-t6/certificate.json','reductions':[{'t':q,'identityVerified':True,'frobeniusDistanceSquared':'0','purity':str(Q(1,q+1)),'evidence':'Certified rho_6=I_7/7 followed by partial trace.'} for q in range(1,7)],'nextOrderWitness':{'order':7,'row':0,'column':0,'nonzeroEntryProven':True,'strictInequalityLatex':r'(\rho_7)_{00}-\frac18<-\frac1{40000}','method':'Exact rational interval around the isolating-root center; the center value is not claimed to equal the root value.'}}}
 v=Q(report['nextOrderWitness']['scaledCenterDefect']);e=Q(report['nextOrderWitness']['scaledRootErrorBound']);assert v+e<-Q(1,40000)
 states.append(n28);states.sort(key=lambda x:x['n'])
 bounds_path=HERE/'proofs.json'
 if bounds_path.exists():
  bounds={row['n']:row for row in json.loads(bounds_path.read_text())['upperBounds']}
  for rec in states:
   bound=bounds.get(rec['n'])
   if bound and rec['achievedOrder']==bound['excludedOrder']-1:
    rec['verification'].update(globalMaximumOrderProven=True,maximumOrderUpperBound=rec['achievedOrder'],maximumOrderProofId=bound['proofId'])
    rec['notes']=f"Exactly order {rec['achievedOrder']}, not order {rec['achievedOrder']+1}. A separate universal exclusion of the next order proves this is the globally maximum pure AC order at N={rec['n']}."
    rec['references'].append({'title':'Universal next-order exclusion proof','url':f"./research/n27-n30/n{rec['n']}-t{bound['excludedOrder']}/proof.tex"})
 (HERE/'states.json').write_text(json.dumps({'states':states},indent=2,ensure_ascii=False)+'\n')
 print('Built',[(x['n'],x['achievedOrder'],len(x['majoranaPoints'])) for x in states])
if __name__=='__main__':main()
