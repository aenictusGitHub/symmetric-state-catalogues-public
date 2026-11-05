#!/usr/bin/env python3
"""Reproduce the exact sparse AC4 witness and its catalogue import records."""
from pathlib import Path
import importlib.util,json
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
ver=load('n19_exact',HERE/'verify.py');pure=load('pure_states',ROOT/'anticoherent-state-catalogue/research/build_states.py')
c={r:s.Add(*(s.Rational(q.numerator,q.denominator)*s.I**e*s.sqrt(rad)for(e,rad),q in v.items()))for r,v in ver.C.items()}
references=[{'title':'Exact sparse AC4 construction and verification','url':'../ac-subspace-catalogue/research/n19-t4/index.html'},{'title':'Exact Dicke coefficients','url':'../ac-subspace-catalogue/research/n19-t4/basis.json'},{'title':'Independent exact reduced-state verification','url':'../ac-subspace-catalogue/research/n19-t4/verification.json'}]
references.append({'title':'Exact maximum AC4 dimension: plane exclusion','url':'../ac-subspace-catalogue/research/n19-t4/exclusion/index.html'})
st=pure.state(19,4,'Five-component analytical fourth-order state at N=19',c,references,provenance='Independent catalogue construction: rational diagonal moment quadrature and an exact phase triangle.',notes='A new exact AC4 witness, not asserted to equal the earlier numerical C3 state. The independent four-copy exclusion proves that no pure AC5 state exists on nineteen qubits, so the maximum pure-state AC order is four. A separate exact Young(2,2) positive-moment certificate excludes every AC4 plane, proving K(19,4)=1.',id='certified-n19-t4')
stateLatex=r'|\psi\rangle=\frac{\sqrt{194}}{40}|D_{19}^{(0)}\rangle+A|D_{19}^{(6)}\rangle+\frac{\sqrt{266}}{160}|D_{19}^{(10)}\rangle+\frac{\sqrt{7410}}{160}|D_{19}^{(14)}\rangle+B|D_{19}^{(18)}\rangle'
definitions=r'\begin{aligned}A&=-\frac{9673\sqrt{17290}}{5605600}-i\frac{\sqrt{2918881141085}}{2802800},\\B&=-\frac{2837\sqrt{323323}}{10210200}-i\frac{\sqrt{218332309353158}}{40840800}.\end{aligned}'
base={'id':'radical-n19-k1-t4-basis-1','n':19,'label':'Exact sparse fourth-order state','stateLatex':stateLatex,'coefficientDefinitionLatex':definitions,'coefficients':st['coefficients'],'exactCoefficients':st['exactCoefficients'],'evidence':'Exact AC4 witness. Majorana coordinates numerically render the exact Dicke coefficients.'}
source={'id':'radical-n19-k1-t4','n':19,'t':4,'dimension':1,'label':'Five-component radical AC4 line','evidence':'exact','references':[{'label':r['title'],'url':r['url'].replace('../ac-subspace-catalogue/','')}for r in references],'notes':st['notes'],'derivationLatex':r'(p_0,p_6,p_{10},p_{14},p_{18})=\left(\frac{97}{800},\frac{1083}{2560},\frac{133}{12800},\frac{741}{2560},\frac{399}{2560}\right)','basis':[base],'verification':{'exact':True,'method':'Every reduced matrix entry through order4 checked exactly with independent standard-library squarefree-radical arithmetic.','orthogonalityError':0,'crossReductionError':0,'fifthOrderFrobeniusDistanceSquared':'8742251/10653696000'}}
(HERE/'basis.json').write_text(json.dumps(source,indent=2)+'\n')
st['stateLatex']=r'\begin{gathered}'+stateLatex+r'\\[4pt]'+definitions+r'\end{gathered}';st['sourceSubspaceId']=source['id'];st['verification']['globalMaximumOrderProven']=True;st['references'].append({'title':'Universal fifth-order exclusion: four-copy proof','url':'../ac-subspace-catalogue/research/n19-t5/index.html'})
(HERE/'pure-import.json').write_text(json.dumps({'states':[st],'geometryImports':{},'method':'Independent exact rational-moment/phase-triangle construction; all normalized-Dicke reductions independently verified.'},indent=2)+'\n')
report=ver.verify();(HERE/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print('Exact witness and import ready',flush=True)
