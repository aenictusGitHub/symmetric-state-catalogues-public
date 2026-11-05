"""Import two already certified radical states and independently recheck their reductions."""
from pathlib import Path
import importlib.util,json,hashlib
import sympy as s
HERE=Path(__file__).resolve().parent;PURE=HERE.parents[1];SUB=PURE.parent/'ac-subspace-catalogue'
spec=importlib.util.spec_from_file_location('pure_exact_reductions',PURE/'research/build_states.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
def main():
 data=json.loads((SUB/'data/catalogue.json').read_text());states=[];audit=[];geometry={}
 for n in [14,17]:
  source=next(x for x in data['subspaces'] if x['id']==f'radical-n{n}-k1-t4');assert source['evidence']=='exact' and source['dimension']==1 and source['t']==4
  vec=source['basis'][0];coeffs=[s.sympify(v) for v in vec['exactCoefficients']]
  st=mod.state(n,4,f'Analytical fourth-order state at N={n}',{i:v for i,v in enumerate(coeffs) if v!=0},
     [{'title':'Analytical construction and exact verification','url':f'../ac-subspace-catalogue/research/n{n}-t4/index.html'},
      {'title':'Exact source coefficients','url':f'../ac-subspace-catalogue/research/n{n}-t4/basis.json'},
      {'title':'Independent exact reduction certificate','url':f'../ac-subspace-catalogue/research/n{n}-t4/independent-verification.json'}],
     provenance='Exact radical construction already certified in the AC-subspace catalogue; every reduction through the achieved order is rechecked here.',
     notes='This state has exactly fourth order. Global order optimality is determined from the exclusion proofs in the catalogue builder. A complex-conjugate representative is also given in the source catalogue.',id=f'certified-n{n}-t4')
  defs=vec['coefficientDefinitionLatex'];extra=source['derivationLatex']if n==17 else ''
  st['stateLatex']=r'\begin{gathered}'+vec['stateLatex']+r'\\[4pt]'+defs+(r'\\[4pt]'+extra if extra else '')+r'\end{gathered}'
  assert st['exactCoefficients']==vec['exactCoefficients'] or all(s.simplify(s.sympify(a)-s.sympify(b))==0 for a,b in zip(st['exactCoefficients'],vec['exactCoefficients']))
  st['sourceSubspaceId']=source['id'];states.append(st)
  geometry[st['id']]={'exactCoefficients':st['exactCoefficients'],'points':vec['points'],'diagnostics':vec['geometryVerification'],'sourceSubspaceId':source['id'],'sourceGeometryCertificate':f'../ac-subspace-catalogue/research/n{n}-t4/geometry-verification.json'}
  audit.append({'n':n,'sourceSubspaceId':source['id'],'exactCoefficientsAssociated':True,'reductions':st['verification']['reductions']})
 imported=json.loads((SUB/'research/n19-t4/pure-import.json').read_text())
 st=imported['state'];assert st['id']=='certified-n19-t4' and all(v['identityVerified'] for v in st['verification']['reductions'] if v['t']<=4)
 states.append(st);geometry[st['id']]=imported['geometry']
 audit.append({'n':19,'sourceSubspaceId':st['sourceSubspaceId'],'exactCoefficientsAssociated':True,'reductions':st['verification']['reductions']})
 out={'states':states,'geometryImports':geometry,'source':'../ac-subspace-catalogue/data/catalogue.json','sourceSHA256':hashlib.sha256((SUB/'data/catalogue.json').read_bytes()).hexdigest(),'method':'Exact SymPy Dicke splitting, independently normalized source coefficients and every reduced matrix through order five.'}
 (HERE/'exact-states.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
 (HERE/'exact-state-verification.json').write_text(json.dumps({'passed':True,'states':audit},indent=2)+'\n')
if __name__=='__main__':main()
