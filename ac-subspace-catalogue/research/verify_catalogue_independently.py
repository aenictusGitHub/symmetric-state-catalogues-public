#!/usr/bin/env python3
"""Independent scientific audit of the saved AC-subspace catalogue.
Requires Python 3, NumPy and SymPy. No source code from the builder is imported.
"""
from pathlib import Path
import json,hashlib,math,importlib.util,os
from urllib.parse import urlsplit,unquote
import numpy as np
import sympy as s
from sympy.physics.wigner import wigner_6j
ROOT=Path(__file__).resolve().parents[1]

def verify_general_entry_bound(n,t):
 """Independent all-shift counting and exact equality tests; no bounds.py import."""
 counts=[sum(l2+1 for l2 in range(h,t+1,2)) for h in range(t+1)]
 upper=min((n-h+1)//counts[h] for h in range(t+1));packing=upper
 odd=t%2==1;r=t//2;d=(r+1)*(r+2) if odd else (r+1)**2
 excluded=False
 if upper and upper*d==(n if odd else n+1):
  j=s.Rational(n,2);jp=j-s.Rational(1,2) if odd else j
  def eigen(L):
   return s.simplify(sum((-1)**(n+L+a)*(2*a+(2 if odd else 1))*wigner_6j(j,j,L,jp,jp,s.Rational(2*a+1,2) if odd else a) for a in range(r+1)))
  roots=[L for L in range(t+1,n if odd else n+1) if eigen(L)==0]
  excluded=not roots and (not odd or n>=7)
  if odd and upper==1 and r>=1 and eigen(t+1)!=0:excluded=True
  if excluded:upper-=1
 return packing,upper,excluded

def verify_algebraic_cases():
 """Recheck the exact certificates, then associate their bases with the export."""
 cases={}
 for n in [14,15,17,19,20,24,35,42]:
  folder=ROOT/f'research/t2-nesting-bounds/n{n}-k{5 if n==35 else 6}' if n in [35,42] else ROOT/'research/t2-general/n24-k5' if n==24 else ROOT/f'research/n{n}-t2' if n in [19,20] else ROOT/f'research/algebraic-k3/n{n}'
  verifier=ROOT/'research/t2-nesting-bounds/verify_parity_independently.py' if n in [35,42] else folder/('verify_independently.py' if n==24 else 'verify.py')
  spec=importlib.util.spec_from_file_location(f'ac_k3_exact_verifier_{n}',verifier)
  module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
  report=module.audit(folder) if n in [35,42] else module.verify();assert report['passed']
  if n in [14,35,42]:count=report['fullCrossReductionEntriesChecked']
  elif n in [19,20,24]:count=report['complexCrossEntriesExactlyVerified']
  elif n==15:count=report['certificate']['fullCrossVerification']['complexEntries']
  else:count=report['completeRealImaginaryReducedMatrixComponents']//2
  dimension=6 if n==42 else 5 if n in [24,35] else 4 if n in [19,20] else 3
  assert count==9*dimension*dimension
  raw=(folder/'basis.json').read_bytes();source=json.loads(raw)
  if 'subspaces' in source:source=source['subspaces'][0]
  assert source['n']==n and source['dimension']==dimension and source.get('t',source.get('order'))==2
  cases[source['id']]={'source':source,'report':{'n':n,'passed':True,'complexCrossEntriesExactlyVerified':count,'basisSHA256':hashlib.sha256(raw).hexdigest(),'method':'Direct exact radical field reduction.' if n==14 else 'Exact rational contraction certificate and all reduced-matrix polynomials reconstructed from Dicke splitting.'}}
 return cases

def verify_four_dimensional_geometry(geometry_module,source):
 """Verify the rational root disks and associate them with the certified basis."""
 n=source['n'];assert n in [19,20,24,35,42]
 k=source['dimension']
 folder=ROOT/f'research/t2-nesting-bounds/n{n}-k{k}' if n in [35,42] else ROOT/'research/t2-general/n24-k5' if n==24 else ROOT/f'research/n{n}-t2'
 spec=importlib.util.spec_from_file_location(f'n{n}_exact_geometry_verifier',folder/'verify_geometry.py')
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 report=module.verify();assert report['passed'] and report['recordCount']==k
 saved=json.loads((folder/'geometry.json').read_text());records=saved['records']
 expected_id=f'algebraic-parity-n{n}-k{k}-t2' if n in [35,42] else f'algebraic-n{n}-k{k}-t2'
 assert len(records)==k and source['id']==expected_id
 parameters=source['algebraicDefinition']['parameterIntervals']
 for index,row in enumerate(source['basis']):
  record=records[source['id']+'-basis-'+str(index+1)]
  assert record['sourceExactCoefficients']==row['exactCoefficients']
  assert record['sourceParameterIntervals']==parameters
  assert record['sourceSHA256']==geometry_module.hashes(row,parameters)
 return records,report

def verify_n12_exclusion():
 """Reconstruct the exact pure-complement identity, independently of the builder."""
 cert=json.loads((ROOT/'research/n12-t2/certificate.json').read_text())
 js=list(range(0,13,2))
 A=s.Matrix(13,7,lambda L,c:(2*L+1)*wigner_6j(6,6,L,6,6,js[c]))
 beta=s.Matrix([[3]+[s.Rational(28,L*(L+1)-112) for L in range(1,13)]])
 f=cert['identity'];rhs=s.Rational(f['constant'])*s.ones(1,7)
 rhs+=sum([s.Rational(f['a'+str(L)])*A[L,:] for L in range(1,4)],s.zeros(1,7))
 rhs+=s.Matrix([list(map(s.Rational,f['pJ']))])
 assert -beta*A==rhs
 assert s.Rational(f['constant'])>0 and s.Rational(f['a3'])>0
 assert all(s.Rational(v)>=0 for v in f['pJ'])
 return {'verifiedExactly':True,'certificateSHA256':hashlib.sha256((ROOT/'research/n12-t2/certificate.json').read_bytes()).hexdigest(),'rankOneComplement':13-4*3,'spinErrorNormalization':14,'evenSpinSectorsChecked':7,'conclusion':'No rank-three second-order AC subspace at N=12.'}

def split(c,n,t):
 return s.Matrix(t+1,n-t+1,lambda a,b:c[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(n-t,b)/s.binomial(n,a+b)))

def symbolic(rows,n,t):
 vectors=[s.Matrix([s.sympify(z) for z in row]) for row in rows];B=s.Matrix.hstack(*vectors)
 gram=B.H*B;checked=0
 for i in range(B.cols):
  for j in range(B.cols):assert s.simplify(s.expand_complex(gram[i,j]-(1 if i==j else 0)))==0,('Gram',i,j)
 matrices=[split(c,n,t) for c in vectors]
 for i,A in enumerate(matrices):
  for j,D in enumerate(matrices):
   deviation=A*D.H-(s.eye(t+1)/(t+1) if i==j else s.zeros(t+1))
   for a in range(t+1):
    for b in range(t+1):
     checked+=1;value=s.simplify(s.expand_complex(deviation[a,b]));assert value==0,('cross',i,j,a,b,str(value))
 return checked

def numerical(rows,n,t):
 B=np.array([[complex(*z) for z in row] for row in rows],complex).T;gram=B.conj().T@B
 mats=[np.array([[B[a+b,i]*math.sqrt(math.comb(t,a)*math.comb(n-t,b)/math.comb(n,a+b)) for b in range(n-t+1)] for a in range(t+1)]) for i in range(B.shape[1])]
 error=max(float(np.linalg.norm(A@D.conj().T-gram[j,i]*np.eye(t+1)/(t+1))) for i,A in enumerate(mats) for j,D in enumerate(mats))
 return {'orthonormalityError':float(np.linalg.norm(gram-np.eye(B.shape[1]))),'maxCrossReductionFrobeniusError':error,'rank':int(np.linalg.matrix_rank(B))}

def main():
 raw=(ROOT/'data/catalogue.json').read_bytes();data=json.loads(raw)
 n12_exclusion=verify_n12_exclusion()
 spec=importlib.util.spec_from_file_location('n9_exact_wedge_exclusion',ROOT/'research/n9-t2/verify.py')
 n9_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n9_module)
 n9_exclusion=n9_module.verify();assert n9_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n11_exact_wedge_exclusion',ROOT/'research/n11-t1/verify.py')
 n11_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n11_module)
 n11_exclusion=n11_module.verify();assert n11_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n13_exact_wedge_exclusion',ROOT/'research/n13-t1/verify.py')
 n13_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n13_module)
 n13_exclusion=n13_module.verify();assert n13_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n13_conditional_exclusion',ROOT/'research/n13-t2/containing-u/verify.py')
 conditional_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(conditional_module)
 n13_conditional=conditional_module.verify();assert n13_conditional['passed']
 spec=importlib.util.spec_from_file_location('third_order_spin5half_inertia',ROOT/'research/t3-spin5half-inertia/verify_independently.py')
 spin5half_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(spin5half_module)
 spin5half_inertia=spin5half_module.verify();assert spin5half_inertia['passed']
 spec=importlib.util.spec_from_file_location('n20_sixth_order_exclusion',ROOT/'research/n20-t6/verify_independently.py')
 n20_sixth_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n20_sixth_module)
 n20_sixth_exclusion=n20_sixth_module.verify();assert n20_sixth_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n18_fourth_order_inertia',ROOT/'research/n18-t4/verify_independently.py')
 n18_fourth_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n18_fourth_module)
 n18_fourth_inertia=n18_fourth_module.verify();assert n18_fourth_inertia['passed']
 spec=importlib.util.spec_from_file_location('n18_fourth_order_witness',ROOT/'research/n18-t4/verify_witness.py')
 n18_fourth_lower_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n18_fourth_lower_module)
 n18_fourth_lower=n18_fourth_lower_module.verify();assert n18_fourth_lower['passed'] and n18_fourth_lower['lowerBound']==1
 spec=importlib.util.spec_from_file_location('n19_fifth_order_exclusion',ROOT/'research/n19-t5/verify_independently.py')
 n19_fifth_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n19_fifth_module)
 n19_fifth_exclusion=n19_fifth_module.verify();assert n19_fifth_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n19_sixth_order_exclusion',ROOT/'research/n19-t6/verify_independently.py')
 n19_sixth_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n19_sixth_module)
 n19_sixth_exclusion=n19_sixth_module.verify();assert n19_sixth_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n18_fifth_order_witness',ROOT/'research/n18-t5/verify.py')
 n18_fifth_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n18_fifth_module)
 n18_fifth_witness=n18_fifth_module.verify();assert n18_fifth_witness['passed'] and n18_fifth_witness['maximumDimension']==1
 spec=importlib.util.spec_from_file_location('n18_sixth_order_exclusion',ROOT/'research/n18-t6/verify_independently.py')
 n18_sixth_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n18_sixth_module)
 n18_sixth_exclusion=n18_sixth_module.verify();assert n18_sixth_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n17_fifth_order_ppt_exclusion',ROOT/'research/n17-t5/verify_independently.py')
 n17_fifth_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n17_fifth_module)
 n17_fifth_exclusion=n17_fifth_module.verify();assert n17_fifth_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n17_sixth_order_complement',ROOT/'research/n17-t6/verify_channel.py')
 n17_sixth_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n17_sixth_module)
 n17_sixth_exclusion=n17_sixth_module.verify();assert n17_sixth_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n17_rank16_exclusion',ROOT/'research/n17-t4/verify_independently.py')
 n17_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n17_module)
 n17_exclusion=n17_module.exclusion();assert n17_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n16_sixth_order_independent_exclusion',ROOT/'research/exact-maxima-patterns/even-tight-independent-verify.py')
 n16_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n16_module)
 n16_exclusion=n16_module.verify();assert n16_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n13_fourth_order_trace_exclusion',ROOT/'research/n13-t4/verify_independently.py')
 n13_fourth_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n13_fourth_module)
 n13_fourth_exclusion=n13_fourth_module.verify();assert n13_fourth_exclusion['passed'] and n13_fourth_exclusion['traceCertificate']['traceAnnihilatesWholeKernelExactly']
 spec=importlib.util.spec_from_file_location('n15_fourth_order_ppt_exclusion',ROOT/'research/n15-t4/verify_independently.py')
 n15_fourth_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n15_fourth_module)
 n15_fourth_exclusion=n15_fourth_module.verify();assert n15_fourth_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n15_third_order_exclusion',ROOT/'research/n15-t3/verify_independently.py')
 n15_third_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n15_third_module)
 n15_third_exclusion=n15_third_module.verify();assert n15_third_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n14_third_order_exclusion',ROOT/'research/n14-t3/verify_independently.py')
 n14_third_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n14_third_module)
 n14_third_exclusion=n14_third_module.verify();assert n14_third_exclusion['passed']
 spec=importlib.util.spec_from_file_location('n14_fifth_order_exclusion',ROOT/'research/n14-t5/verify_independently.py')
 n14_fifth_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(n14_fifth_module)
 n14_fifth_exclusion=n14_fifth_module.verify();assert n14_fifth_exclusion['passed'] and n14_fifth_exclusion['positiveDiagonalCertificate']['formsHaveRankTwo']
 n8_proofs={}
 for order in [1,2]:
  spec=importlib.util.spec_from_file_location(f'n8_exact_order{order}',ROOT/f'research/n8/t{order}/verify.py')
  module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
  n8_proofs[order]=module.verify();assert n8_proofs[order]['passed']
 plane_cases={};fourth_exclusions={}
 for n in [19,20]:
  folder=ROOT/f'research/n{n}-t3'
  spec=importlib.util.spec_from_file_location(f'n{n}_third_order_plane',folder/'verify_independently.py')
  module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
  report=module.verify();assert report['passed'] and report['fullCrossReductionEntriesChecked']==64
  source=json.loads((folder/'basis.json').read_text())
  plane_cases[source['id']]={'source':source,'report':report}
  spec=importlib.util.spec_from_file_location(f'n{n}_fourth_order_plane_exclusion',ROOT/f'research/n{n}-t4/exclusion/verify_independently.py')
  module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
  if os.environ.get('AC_REUSE_FOURTH_ORDER_AUDITS')=='1':
   report=json.loads((ROOT/f'research/n{n}-t4/exclusion/independent-verification.json').read_text())
   digest=hashlib.sha256((ROOT/f'research/n{n}-t4/exclusion/certificate.json').read_bytes()).hexdigest()
   expected={19:'9c3ab66655eb78b3e5a723164f2cf26793830543d593c713d9623f330349558d',20:'543b47263c8deb5d30677bcab8793a305426853bbbb433836cc3d95cde6b052c'}
   assert digest==expected[n]==report['certificateSHA256']
   assert report['physicalConstraintRowsIndependentlyReconstructed'] and report['noFloatingPointProofSteps']
   assert report['constraintRowsReconstructed']==(3535 if n==19 else 4253)
   assert report['positiveDefiniteRowCertificate']['allLDLReconstructionsExact']
  else:report=module.verify()
  assert report['passed'];fourth_exclusions[str(n)]=report
 algebraic_cases=verify_algebraic_cases()
 algebraic_cases.update(plane_cases)
 geometry_spec=importlib.util.spec_from_file_location('ac_k3_geometry_verifier',ROOT/'research/algebraic-k3/verify_geometry.py')
 geometry_module=importlib.util.module_from_spec(geometry_spec);geometry_spec.loader.exec_module(geometry_module)
 geometry_report=geometry_module.verify();assert geometry_report['passed']
 algebraic_geometry=json.loads((ROOT/'research/algebraic-k3/geometry.json').read_text())['records']
 plane_geometry_reports={}
 for sid,case in plane_cases.items():
  n=case['source']['n'];folder=ROOT/f'research/n{n}-t3'
  spec=importlib.util.spec_from_file_location(f'n{n}_third_order_geometry',folder/'verify_geometry.py')
  module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
  report=module.verify();assert report['passed'] and report['recordCount']==2
  records=json.loads((folder/'geometry.json').read_text())['records']
  for row in case['source']['basis']:assert records[row['id']]['sourceExactCoefficients']==row['exactCoefficients']
  assert not set(records)&set(algebraic_geometry);algebraic_geometry.update(records);plane_geometry_reports[str(n)]=report
 n20_source=algebraic_cases['algebraic-n20-k4-t2']['source']
 n20_geometry,n20_geometry_report=verify_four_dimensional_geometry(geometry_module,n20_source)
 assert not set(n20_geometry)&set(algebraic_geometry)
 algebraic_geometry.update(n20_geometry)
 n19_geometry,n19_geometry_report=verify_four_dimensional_geometry(geometry_module,algebraic_cases['algebraic-n19-k4-t2']['source'])
 assert not set(n19_geometry)&set(algebraic_geometry)
 algebraic_geometry.update(n19_geometry)
 n24_geometry,n24_geometry_report=verify_four_dimensional_geometry(geometry_module,algebraic_cases['algebraic-n24-k5-t2']['source'])
 assert not set(n24_geometry)&set(algebraic_geometry)
 algebraic_geometry.update(n24_geometry)
 parity_geometry_reports={}
 for n,k in [(35,5),(42,6)]:
  source=algebraic_cases[f'algebraic-parity-n{n}-k{k}-t2']['source']
  records,parity_report=verify_four_dimensional_geometry(geometry_module,source)
  assert not set(records)&set(algebraic_geometry)
  algebraic_geometry.update(records);parity_geometry_reports[str(n)]=parity_report
 pure=json.loads((ROOT/'research/pure-catalogue.json').read_text());pure_entries={(e['n'],e['t']):e for e in pure['entries']}
 paper=json.loads((ROOT/'research/paper/findings.json').read_text());paper_by_id={e['id']:e for e in paper['exactWitnesses']}
 proofids={p['id'] for p in data['proofs']};byid={r['id']:r for r in data['subspaces']};assert len(byid)==len(data['subspaces'])
 exact=[];numeric=[];scalar_count=0;max_exact_float_error=0;basis_count=0
 for rec in data['subspaces']:
  n,t,k=rec['n'],rec['t'],rec['dimension'];assert k==len(rec['basis']) and k>=1 and n>=t>=1
  assert rec['evidence'] in ['exact','numerical'] and rec['verification']['exact']==(rec['evidence']=='exact')
  rows=[b['coefficients'] for b in rec['basis']]
  for b in rec['basis']:
   basis_count+=1;assert len(b['coefficients'])==n+1 and b['n']==n
   assert b['points'] and sum(p['multiplicity'] for p in b['points'])==n
   assert all(p['multiplicity']>=1 and isinstance(p['multiplicity'],int) for p in b['points'])
   if rec['id'] in algebraic_cases:
    assert b['points']==algebraic_geometry[b['id']]['points']
    assert b['geometryVerification']==algebraic_geometry[b['id']]['diagnostics']
  check=numerical(rows,n,t);assert check['rank']==k
  if rec['evidence']=='exact':
   assert all(len(b.get('exactCoefficients',[]))==n+1 for b in rec['basis'])
   expressions=[b['exactCoefficients'] for b in rec['basis']]
   parameters={};center={}
   if rec['id'] in algebraic_cases:
    certified=algebraic_cases[rec['id']];source=certified['source']
    assert expressions==[b['exactCoefficients'] for b in source['basis']]
    if source.get('exactnessKind')=='isolated-algebraic-root':
     assert rec['algebraicDefinition']==source['algebraicDefinition']
     intervals=source['algebraicDefinition']['parameterIntervals']
     parameters={name:s.Symbol(name,real=True) for name in intervals}
     center={parameters[name]:(s.Rational(bounds[0])+s.Rational(bounds[1]))/2 for name,bounds in intervals.items()}
    count=certified['report']['complexCrossEntriesExactlyVerified']
   else:count=symbolic(expressions,n,t)
   scalar_count+=count
   for values,coeffs in zip(expressions,rows):
    err=max(abs(complex(s.N(s.sympify(x,locals=parameters).subs(center),25))-complex(*z)) for x,z in zip(values,coeffs));assert err<3e-12;max_exact_float_error=max(max_exact_float_error,err)
   exact.append({'id':rec['id'],'n':n,'t':t,'dimension':k,'orthonormal':True,'allOrderedCrossEntriesExactlyVerified':True,'scalarEntriesChecked':count})
  else:
   assert not any('exactCoefficients' in b for b in rec['basis'])
   assert check['orthonormalityError']<1e-7 and check['maxCrossReductionFrobeniusError']<1e-7,(rec['id'],check)
   numeric.append({'id':rec['id'],'n':n,'t':t,'dimension':k,**check,'exactExistenceClaim':False})
  if rec['id'] in paper_by_id:
   source=paper_by_id[rec['id']];assert n==source['n'] and k==source['dimension'] and t==source.get('verifiedAchievedOrder',source['order'])
   assert expressions==source['exactBasis']
  print(rec['id'],'verified',flush=True)
 assert set(paper_by_id)<=set(byid),'missing paper witnesses'
 entrychecks=[];zero_exclusions=0;positive_tight=0;conjecture_count=0;saturation_exclusions=[]
 for entry in data['entries']:
  n,t=entry['n'],entry['t'];assert n in data['numbers'] and t in data['orders'] and t<=n
  eligible=[r for r in data['subspaces'] if r['n']==n and r['t']>=t];rigorous=[r for r in eligible if r['evidence']=='exact'];numerics=[r for r in eligible if r['evidence']=='numerical']
  lower=max((r['dimension'] for r in rigorous),default=0)
  schmidt=(n-t+1)//(t+1);gram=(n+1)//((t//2+1)**2);upper=max(0,min(schmidt,gram))
  if n>=2:
   # Independent spin-one coupling calculation.  Since t>=1, each
   # compressed J component vanishes on the whole code; its tensor
   # product with the three-dimensional ancillary space is totally
   # isotropic for J.S.  Count both strict spectral signs.
   spin=s.Rational(n,2)
   total_spins=[spin-1,spin,spin+1]
   eigenvalues=[s.factor((L*(L+1)-spin*(spin+1)-2)/2) for L in total_spins]
   positive=sum(int(2*L+1) for L,value in zip(total_spins,eigenvalues) if value>0)
   negative=sum(int(2*L+1) for L,value in zip(total_spins,eigenvalues) if value<0)
   assert all(eigenvalues) and (positive,negative)==(n+3,2*n)
   inertia_upper=min(positive,negative)//3
   assert entry['firstOrderInertiaUpper']==inertia_upper
   if t==1:assert 'first-order-optimality' in entry['proofIds']
   upper=min(upper,inertia_upper)
  else:
   # At N=1 the j-1 summand is absent.  J.S has eigenvalues 1/2
   # and -1 with multiplicities 4 and 2, hence 3k<=2 and k=0.
   assert n==1 and entry['firstOrderInertiaUpper']==0
   assert 'first-order-optimality' in entry['proofIds']
   upper=0
  if t==1 and n%2==0 and n>=8:
   assert n8_proofs[1]['passed'] and 'first-order-saturation-obstruction' in entry['proofIds']
   if n>=10:assert (n//2)**2>2*n+1
   upper=min(upper,n//2-1)
  if n==8 and t>=2:
   assert n8_proofs[2]['passed'] and 'n8-second-order-complement' in entry['proofIds']
   upper=min(upper,1)
  if n==9 and t==2:
   assert n9_exclusion['passed'] and 'n9-second-order-wedge' in entry['proofIds']
   upper=min(upper,1)
  if n==11 and t==1:
   assert n11_exclusion['passed'] and 'n11-first-order-wedge' in entry['proofIds']
   upper=min(upper,4)
  if n==13 and t==1:
   assert n13_exclusion['passed'] and 'n13-first-order-wedge' in entry['proofIds']
   upper=min(upper,5)
  general_packing,general_upper,general_excluded=verify_general_entry_bound(n,t)
  assert 'rectangular-tensor-bound' in entry['proofIds']
  if general_excluded:assert 'tensor-packing-saturation' in entry['proofIds']
  assert entry['generalBound']['packingUpper']==general_packing
  assert entry['generalBound']['upper']==general_upper
  assert entry['generalBound']['saturationExcluded']==general_excluded
  upper=min(upper,general_upper)
  if n>=8 and t>=2:
   j=s.Rational(n,2);aux=s.Rational(3,2)
   spins=[j-aux+i for i in range(4)]
   xs=[(L*(L+1)-j*(j+1)-aux*(aux+1))/2 for L in spins]
   eig=[s.factor(x*x-s.Rational(n-4,2)*x-s.Rational(5*n*(n+2),16)) for x in xs]
   assert eig[0]>0 and all(v<0 for v in eig[1:])
   positive=sum(int(2*L+1) for L,v in zip(spins,eig) if v>0)
   assert positive==n-2 and entry['ancillaInertiaUpper']==positive//4
   assert 'ancilla-inertia-bound' in entry['proofIds']
   upper=min(upper,positive//4)
  else:assert entry['ancillaInertiaUpper'] is None
  if n>=4 and t>=2:
   # Independent spin-two spectrum of X=J.S and its quadratic witness.
   # Whole-subspace 2-AC implies compression(X)=0 and
   # compression(X^2)=2*j*(j+1) I, so the five-fold code is isotropic.
   j=s.Rational(n,2)
   spins=[j-2+a for a in range(5)]
   xs=[(L*(L+1)-j*(j+1)-6)/2 for L in spins]
   eig=[s.factor(x*x+(j+3)*x-2*j*(j+1)) for x in xs]
   assert all(v<0 for v in eig[:4]) and eig[4]>0
   positive=sum(int(2*L+1) for L,v in zip(spins,eig) if v>0)
   negative=sum(int(2*L+1) for L,v in zip(spins,eig) if v<0)
   assert (positive,negative)==(n+5,4*n)
   spin2_upper=min(positive,negative)//5
   assert spin2_upper==(n+5)//5 and entry['spin2InertiaUpper']==spin2_upper
   assert 'second-order-spin2-inertia' in entry['proofIds']
   upper=min(upper,spin2_upper)
  else:assert entry['spin2InertiaUpper'] is None
  if n>=6 and t>=3:
   # Reconstruct the cubic spectrum and its normalized trace independently.
   j=s.Rational(n,2);spins=[j-3+a for a in range(7)]
   xs=[(L*(L+1)-j*(j+1)-12)/2 for L in spins]
   values=[s.factor((x+2*j+5)*(x+j+6)*(x-2*j+3)-15*(2*j*j+3*j+6)) for x in xs]
   assert sum((2*L+1)*v for L,v in zip(spins,values))==0
   assert all(v<0 for v in values[:-1]) and values[-1]>0
   assert int(2*spins[-1]+1)==n+7
   assert entry['spin3InertiaUpper']==(n+7)//7 and 'third-order-spin3-inertia' in entry['proofIds']
   upper=min(upper,(n+7)//7)
  else:assert entry['spin3InertiaUpper'] is None
  if n>=9 and t>=4:
   j=s.Rational(n,2);aux=s.Rational(9,2)
   spins=[j-aux+a for a in range(10)]
   xs=[(L*(L+1)-j*(j+1)-aux*(aux+1))/2 for L in spins]
   roots=[xs[i] for i in [0,3,4,8]]
   polynomial=[s.prod(x-root for root in roots) for x in xs]
   mean=s.factor(sum((2*L+1)*v for L,v in zip(spins,polynomial))/(10*(n+1)))
   assert mean>0
   values=[s.factor(v-mean) for v in polynomial]
   assert sum((2*L+1)*v for L,v in zip(spins,values))==0
   assert all(v<0 for v in values[:-1]) and values[-1]>0
   assert int(2*spins[-1]+1)==n+10
   assert entry['spin9HalfInertiaUpper']==(n+10)//10 and 'fourth-order-spin9half-inertia' in entry['proofIds']
   upper=min(upper,(n+10)//10)
  else:assert entry['spin9HalfInertiaUpper'] is None
  if n>=19 and t>=3:
   j=s.Rational(n,2);aux=s.Rational(5,2);spins=[j-aux+a for a in range(6)]
   xs=[(L*(L+1)-j*(j+1)-aux*(aux+1))/2 for L in spins]
   vals=[s.factor(-(2*x+3*j+8)*(2*x-3*j+5)*(2*x-5*j))for x in xs]
   mean=s.factor(sum((2*L+1)*v for L,v in zip(spins,vals))/(6*(n+1)))
   assert mean==s.Rational(20,3)*j*(j-9)*(2*j-1)>0
   eig=[v-mean for v in vals];assert eig[0]>0 and all(v<0 for v in eig[1:])
   assert sum((2*L+1)*v for L,v in zip(spins,eig))==0 and int(2*spins[0]+1)==n-4
   assert entry['thirdOrderSpin5HalfInertiaUpper']==(n-4)//6 and 'third-order-spin5half-inertia' in entry['proofIds']
   upper=min(upper,(n-4)//6)
  else:assert entry['thirdOrderSpin5HalfInertiaUpper'] is None
  if n>=16 and t>=4:
   j=s.Rational(n,2);spins=[j-3+a for a in range(7)]
   xs=[(L*(L+1)-j*(j+1)-12)/2 for L in spins]
   values=[s.factor((x+2*j+5)*(x+6)*(x-j+5)*(x-3*j)) for x in xs]
   mean=s.factor(sum((2*L+1)*v for L,v in zip(spins,values))/(7*(n+1)))
   h=j-8;assert mean==2*j*(4*h**3+60*h*h+239*h+45)>0
   eig=[v-mean for v in values];assert eig[0]>0 and all(v<0 for v in eig[1:])
   assert sum((2*L+1)*v for L,v in zip(spins,eig))==0
   assert int(2*spins[0]+1)==n-5
   assert entry['fourthOrderSpin3InertiaUpper']==(n-5)//7 and 'fourth-order-spin3-inertia' in entry['proofIds']
   upper=min(upper,(n-5)//7)
  else:assert entry['fourthOrderSpin3InertiaUpper'] is None
  if n==18 and t==4:
   assert n18_fourth_inertia['passed'] and n18_fourth_lower['passed'] and 'n18-fourth-order-maximum' in entry['proofIds']
   assert (entry['lower'],entry['upper'],entry['status'])==(1,1,'proven')
  if n==20 and t>=6:
   assert n20_sixth_exclusion['passed'] and 'n20-sixth-order-exclusion' in entry['proofIds']
   upper=0
  if n==19 and t>=5:
   assert n19_fifth_exclusion['passed'] and 'n19-fifth-order-exclusion' in entry['proofIds']
   upper=0
  if n==19 and t>=6:
   assert n19_sixth_exclusion['passed'] and 'n19-sixth-order-exclusion' in entry['proofIds']
   upper=0
  if n==18 and t==5:
   assert n18_fifth_witness['passed'] and 'n18-fifth-order-maximum' in entry['proofIds']
  if n==18 and t>=6:
   assert n18_sixth_exclusion['passed'] and 'n18-sixth-order-exclusion' in entry['proofIds']
   upper=0
  if n==17 and t>=5:
   assert n17_fifth_exclusion['passed'] and 'n17-fifth-order-exclusion' in entry['proofIds']
   upper=0
  if n==17 and t>=6:
   assert n17_sixth_exclusion['passed'] and 'n17-sixth-order-complement' in entry['proofIds']
  if n==16 and t>=6:
   assert n16_exclusion['passed'] and 'n16-sixth-order-exclusion' in entry['proofIds']
   upper=0
  if n==13 and t>=4:
   assert n13_fourth_exclusion['passed'] and 'n13-fourth-order-exclusion' in entry['proofIds']
   upper=0
  if n==15 and t>=3:
   assert n15_third_exclusion['passed'] and 'n15-third-order-maximum' in entry['proofIds']
   upper=min(upper,1)
  if n==15 and t>=4:
   assert n15_fourth_exclusion['passed'] and 'n15-fourth-order-exclusion' in entry['proofIds']
   upper=0
  if n==15 and t==2:
   assert 'n15-second-order-saturation' in entry['proofIds']
   assert not any(L*(L+1)==170 for L in range(1,16))
  if n==14 and t>=3:
   assert n14_third_exclusion['passed'] and 'n14-third-order-maximum' in entry['proofIds']
   upper=min(upper,1)
  if n==14 and t>=5:
   assert n14_fifth_exclusion['passed'] and 'n14-fifth-order-exclusion' in entry['proofIds']
   upper=0
  # Saturation of the four spin-error images requires a surviving tensor rank.
  # The Casimir-map coefficient is 4 - 6 L(L+1)/(n(n+2)).
  if t>=2 and 4*upper==n+1:
   ranks=[L for L in range(1,n+1) if s.Rational(6*L*(L+1),n*(n+2))==4]
   if not ranks:
    upper-=1;assert 'gram-saturation-obstruction' in entry['proofIds']
    saturation_exclusions.append({'n':n,'t':t,'excludedDimension':(n+1)//4,'newUpper':upper,'allowedTensorRanks':ranks})
  if n==12 and t==2:
   assert n12_exclusion['verifiedExactly'] and 'n12-t2-complement' in entry['proofIds'];upper=2
  if n in [19,20] and t>=4:
   assert fourth_exclusions[str(n)]['passed'] and f'n{n}-fourth-order-maximum' in entry['proofIds']
   upper=min(upper,1)
  if n==17 and t>=4:
   assert n17_exclusion['passed'] and 'n17-fourth-order-maximum' in entry['proofIds']
   upper=min(upper,1)
  purecell=pure_entries.get((n,t));pureproof=None
  if purecell and purecell['status']=='impossible':upper=0;pureproof='pure-'+purecell['proofId'];assert pureproof in entry['proofIds'];zero_exclusions+=1
  assert entry['lower']==lower and entry['upper']==upper and lower<=upper,(entry['id'],lower,upper,entry['lower'],entry['upper'])
  assert set(entry['subspaceIds'])=={r['id'] for r in eligible},entry['id']
  assert all(pid in proofids for pid in entry['proofIds'])
  assert {'cross-condition','paper-schmidt-bound','tensor-gram-bound'}<=set(entry['proofIds'])
  expected_status='impossible' if upper==0 else 'proven' if lower==upper else 'range';assert entry['status']==expected_status
  reported=max((r['dimension'] for r in numerics),default=0)
  if n in paper['figure2']['n'] and t in [1,2]:reported=max(reported,paper['figure2']['t'+str(t)][paper['figure2']['n'].index(n)])
  assert (entry['reportedDimension'] or 0)==reported and reported<=upper
  proposed=paper['figure2']['t'+str(t)][paper['figure2']['n'].index(n)] if n in paper['figure2']['n'] and t in [1,2] else None
  assert entry['conjecturedDimension']==proposed
  if proposed is not None:
   conjecture_count+=1
   status='superseded' if lower>proposed else 'confirmed' if lower==upper==proposed else 'attained' if lower==proposed else 'open'
   assert entry['conjectureStatus']==status and entry['conjectureReferences']
   if t==1:assert lower>=proposed,('first-order value not attained',n,proposed,lower)
  else:assert entry['conjectureStatus'] is None and not entry['conjectureReferences']
  if expected_status=='proven':positive_tight+=1
  entrychecks.append({'id':entry['id'],'lower':lower,'upper':upper,'status':expected_status,'numericEvidenceExcludedFromLowerBound':True,'schmidtUpper':schmidt,'tensorGramUpper':gram,'pureExclusion':pureproof})
 assert len({(e['n'],e['t']) for e in data['entries']})==len(data['entries'])
 assert {(e['n'],e['t']) for e in data['entries']}=={(n,t) for n in data['numbers'] for t in data['orders'] if t<=n}
 local_reference_paths=set()
 for group in [data['subspaces'],data['entries'],data['proofs']]:
  for obj in group:
   for ref in obj.get('references',[])+obj.get('reportedReferences',[])+obj.get('conjectureReferences',[]):
    parsed=urlsplit(ref.get('url',''))
    if parsed.scheme or not parsed.path:continue
    target=(ROOT/unquote(parsed.path)).resolve();assert target.exists(),(obj['id'],ref['url']);local_reference_paths.add(str(target))
 known={(19,3):(2,2),(20,3):(2,2),(19,4):(1,1),(20,4):(1,1),(20,6):(0,0),(13,4):(0,0),(13,5):(0,0),(16,6):(0,0),(1,1):(0,0),(7,1):(3,3),(8,1):(3,3),(8,2):(1,1),(8,3):(1,1),(10,1):(4,4),(12,1):(5,5),(10,2):(2,2),(11,2):(2,2),(13,3):(2,2),(12,2):(2,2),(14,2):(3,3),(15,2):(3,3),(17,2):(3,3),(16,2):(3,3),(13,2):(2,2),(18,3):(2,2),(9,3):(0,0),(10,3):(1,1),(11,3):(1,1),(12,3):(1,1),(18,5):(1,1),(20,5):(1,1),(24,5):(1,1),(20,7):(0,0)}
 for key,interval in known.items():
  cell=next(e for e in data['entries'] if (e['n'],e['t'])==key);assert (cell['lower'],cell['upper'])==interval
 n20=next(e for e in data['entries'] if (e['n'],e['t'])==(20,2))
 assert (n20['lower'],n20['upper'],n20['status'])==(4,4,'proven')
 assert 'algebraic-n20-k4-t2' in n20['subspaceIds']
 for n,lower,upper in [(22,3,5),(24,5,5),(35,5,8),(42,6,9)]:
  cell=next(e for e in data['entries'] if (e['n'],e['t'])==(n,2))
  assert (cell['lower'],cell['upper'])==(lower,upper)
  assert 'general-second-order-constructions' in cell['proofIds']
  if n==24:assert 'n24-second-order-maximum' in cell['proofIds'] and cell['status']=='proven'
  if n in [35,42]:assert f'algebraic-parity-n{n}-k{lower}-t2' in cell['subspaceIds'] and cell['status']=='range' and 'nested-parity-second-order-constructions' in cell['proofIds']
 n14=next(e for e in data['entries'] if (e['n'],e['t'])==(14,4))
 assert (n14['lower'],n14['upper'],n14['status'])==(1,1,'proven')
 assert {'radical-n14-k1-t4','line-acliege-n14-t4'}<=set(n14['subspaceIds'])
 assert 'n14-fourth-order-maximum' in n14['proofIds']
 for n,t,k,pid,sids in [
   (16,3,2,'n16-third-order-maximum',{'radical-n16-k2-t3','acliege-Π_16_2_3'}),
   (17,3,2,'n17-third-order-maximum',{'radical-n17-k2-t3','acliege-Π_17_2_3'}),
   (18,2,4,'n18-second-order-maximum',{'radical-n18-k4-t2','acliege-Π_18_4_2'}),
   (19,2,4,'n19-second-order-maximum',{'algebraic-n19-k4-t2','acliege-Π_19_4_2'}),
   (17,4,1,'n17-fourth-order-maximum',{'radical-n17-k1-t4','radical-n17-k1-t4-conjugate','line-acliege-n17-t4'})]:
  cell=next(e for e in data['entries'] if (e['n'],e['t'])==(n,t))
  assert (cell['lower'],cell['upper'],cell['status'])==(k,k,'proven')
  assert sids<=set(cell['subspaceIds']) and pid in cell['proofIds']
 for cell in data['entries']:
  if cell['t']!=1:continue
  n=cell['n'];expected=0 if n==1 else 1 if n in [2,3] else n//3+1
  assert (cell['lower'],cell['upper'])==(expected,expected),('global first-order maximum',n,cell)
 cell=next(e for e in data['entries'] if (e['n'],e['t'])==(9,2))
 assert (cell['lower'],cell['upper'],cell['conjectureStatus'])==(1,1,'confirmed')
 cell=next(e for e in data['entries'] if (e['n'],e['t'])==(13,2))
 assert (cell['lower'],cell['upper'],cell['conjectureStatus'])==(2,2,'confirmed')
 assert 'n13-second-order-inertia' in cell['proofIds']
 assert 'n13-phase-family' in cell['proofIds']
 assert 'n13-sparse-state-exclusion' in cell['proofIds']
 cell=next(e for e in data['entries'] if (e['n'],e['t'])==(13,1))
 assert (cell['lower'],cell['upper'],cell['conjectureStatus'])==(5,5,'confirmed')
 for id in ['family51-n18','sparse-third-n18']:assert byid[id]['t']==4
 report={'passed':True,'catalogueSha256':hashlib.sha256(raw).hexdigest(),'exactSubspacesChecked':len(exact),'numericalSubspacesChecked':len(numeric),'individualBasisVectorsChecked':basis_count,'entriesChecked':len(entrychecks),'paperWitnessesAssociated':len(paper_by_id),'uniqueLocalReferencePathsChecked':len(local_reference_paths),'exactCrossMatrixElementsChecked':scalar_count,'maximumExactFormulaVsSavedCoefficientError':max_exact_float_error,'pureExclusionsApplied':zero_exclusions,'provedPositiveMaxima':positive_tight,'method':'Independent exact normalized-Dicke splitting and full cross-reduction matrices; no builder functions imported. Every lower bound recomputed from exact subspaces only, each upper bound recomputed from Schmidt and tensor-Gram inequalities plus supplied pure-state exclusions, and all proof associations checked.','exactSubspaces':exact,'numericalSubspaces':numeric,'entries':entrychecks,'limitations':['Numerical basis checks remain approximate and do not establish exact existence or maximality.','The linked pure-state exclusion proofs are reused rather than reproved by this audit.','Individual Majorana coordinates are tested separately by reconstructing the vectors.']}
 report['n12RankOneComplementExclusion']=n12_exclusion
 report['n9FourCopyWedgeExclusion']=n9_exclusion
 report['n11FirstOrderExteriorExclusion']=n11_exclusion
 report['n13FirstOrderExteriorExclusion']=n13_exclusion
 report['n13SecondOrderConditionalExclusion']=n13_conditional
 report['n17FourthOrderSaturationLiftExclusion']=n17_exclusion
 report['thirdOrderSpin5HalfInertia']=spin5half_inertia
 report['n19n20ThirdOrderPlanes']={sid:case['report'] for sid,case in plane_cases.items()}
 report['n19n20ThirdOrderPlaneGeometry']=plane_geometry_reports
 report['n19n20FourthOrderExclusions']=fourth_exclusions
 report['n20SixthOrderMatrixExclusion']=n20_sixth_exclusion
 report['n18FourthOrderInertia']=n18_fourth_inertia
 report['n18FourthOrderWitness']=n18_fourth_lower
 report['fourthOrderSpinThreeGeneralBound']={'formula':'7k<=N-5','minimumN':16,'minimumOrder':4,'allApplicableEntriesRecomputed':True}
 report['n19FifthOrderMatrixExclusion']=n19_fifth_exclusion
 report['n19SixthOrderProbabilityExclusion']=n19_sixth_exclusion
 report['n18FifthOrderWitnessAndDimension']=n18_fifth_witness
 report['n18SixthOrderProbabilityExclusion']=n18_sixth_exclusion
 report['n17FifthOrderPPTExclusion']=n17_fifth_exclusion
 report['n17SixthOrderComplementExclusion']=n17_sixth_exclusion
 report['n16SixthOrderExclusion']=n16_exclusion
 report['n13FourthOrderTraceExclusion']=n13_fourth_exclusion
 report['n14FifthOrderPositivityExclusion']=n14_fifth_exclusion
 report['n14ThirdOrderPositivityExclusion']=n14_third_exclusion
 report['n15ThirdOrderPositivityExclusion']=n15_third_exclusion
 report['n15FourthOrderPPTExclusion']=n15_fourth_exclusion
 report['thirdOrderSpinThreeInertia']={'allApplicableEntriesChecked':True,'spectrumAndNormalizedTraceRecomputedIndependently':True,'generalBound':'7k <= N+7 for N>=6 and t>=3'}
 report['fourthOrderSpinNineHalfInertia']={'allApplicableEntriesChecked':True,'spectrumAndNormalizedTraceRecomputedIndependently':True,'generalBound':'10k <= N+10 for N>=9 and t>=4'}
 report['auxiliarySpinInertia']={'allApplicableEntriesChecked':True,'spectrumRecomputedIndependently':True,'generalBound':'4k <= N-2 for N>=8 and t>=2'}
 report['secondOrderSpinTwoInertia']={'allApplicableEntriesChecked':True,'spectrumRecomputedIndependently':True,'generalBound':'5k <= N+5 for N>=4 and t>=2','entriesChecked':sum(e['n']>=4 and e['t']>=2 for e in data['entries'])}
 report['firstOrderSpinOneInertia']={'allApplicableEntriesChecked':True,'spectrumRecomputedIndependently':True,'generalBound':'3k <= min(N+3, 2N) for N>=2 and t>=1','allFirstOrderMaximaExact':True,'firstOrderEntriesChecked':sum(e['t']==1 for e in data['entries'])}
 report['generalRectangularBounds']={'allEntriesRecomputed':True,'method':'All rectangular shifts counted independently; equality multipliers reconstructed with exact six-j coefficients.'}
 report['n8SharpDimensionProofs']=n8_proofs
 assert conjecture_count==33
 report['paperConjecturesChecked']=conjecture_count
 report['allPaperFirstOrderDimensionsAttainedExactly']=True
 report['gramSaturationExclusions']=saturation_exclusions
 report['newThreeDimensionalExactCertificates']=[case['report'] for case in algebraic_cases.values() if case['source']['dimension']==3]
 report['n20FourDimensionalExactCertificate']=algebraic_cases['algebraic-n20-k4-t2']['report']
 report['n20AlgebraicMajoranaVerification']=n20_geometry_report
 report['n19FourDimensionalExactCertificate']=algebraic_cases['algebraic-n19-k4-t2']['report']
 report['n19AlgebraicMajoranaVerification']=n19_geometry_report
 report['n24FiveDimensionalExactCertificate']=algebraic_cases['algebraic-n24-k5-t2']['report']
 report['n24AlgebraicMajoranaVerification']=n24_geometry_report
 report['nestedParityExactCertificates']=[algebraic_cases[f'algebraic-parity-n{n}-k{k}-t2']['report']for n,k in [(35,5),(42,6)]]
 report['nestedParityMajoranaVerification']=parity_geometry_reports
 report['newAlgebraicMajoranaVerification']=geometry_report
 report['allPaperProposedDimensionsAttainedExactly']=all(e['lower']>=e['conjecturedDimension'] for e in data['entries'] if e['conjecturedDimension'] is not None)
 assert report['allPaperProposedDimensionsAttainedExactly']
 report['method']+=' The N=12,t=2 complement certificate is independently reconstructed with exact rational arithmetic.'
 report['method']+=' Spin-error saturation is also excluded whenever no integer tensor rank has the required Casimir eigenvalue. Every conjectured dimension is checked against Figure 2.'
 report['method']+=' The three new dimension-three examples are checked by their exact radical or rational-isolation verifiers; the exported symbolic coefficients and parameter definitions must match the certified basis exactly.'
 report['method']+=' The N=8 first- and second-order exclusions are rechecked exactly, including the auxiliary pure-state recoupling identities and the vectorization lift; its general even-N first-order consequence is applied only on the proved domain.'
 report['method']+=' The N=9 second-order exclusion is reconstructed with exact Young-symmetry contraction equations and a strictly positive diagonal identity.'
 (ROOT/'research/catalogue-independent-audit.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k not in ['exactSubspaces','numericalSubspaces','entries']},indent=2))
if __name__=='__main__':main()
