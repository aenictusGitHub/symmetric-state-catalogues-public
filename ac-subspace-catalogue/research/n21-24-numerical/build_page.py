"""Render the numerical survey without promoting finite searches to proofs."""
from pathlib import Path
from html import escape
import json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def m(t):return '<span data-tex="'+escape(t,quote=True)+'"></span>'
def sci(x):return f'{float(x):.3e}'
def link(url,label):return '<a href="'+escape(url,quote=True)+'">'+escape(label)+'</a>'
def main():
 survey=json.loads((HERE/'bases.json').read_text());audit=json.loads((HERE/'independent-verification.json').read_text());data=json.loads((ROOT/'data/catalogue.json').read_text());entries={(e['n'],e['t']):e for e in data['entries']};checks={e['id']:e for e in audit['candidates']};assert audit['passed']
 b=link('../../index.html#n=21&t=2','← AC-subspace catalogue')+'<h1>Numerical dimensions at N = 21–24</h1><p>High-accuracy searches for whole anticoherent subspaces, with independent checks of every ordered reduced matrix. Failed searches do not prove nonexistence; small residuals do not supply exact existence certificates.</p>'
 b+='<section><h2>Largest dimensions currently found</h2><div class="scroll"><table><thead><tr><th>'+m('N')+'</th>'+''.join('<th>'+m('t='+str(t))+'</th>'for t in range(1,8))+'</tr></thead><tbody>'
 for n in range(21,25):
  b+='<tr><th>'+str(n)+'</th>'
  for t in range(1,8):
   e=entries[n,t];value=max(e['lower'],e['reportedDimension'] or 0);status=e['status']
   text=str(value)if status=='proven' else str(value)+' (num.)'if value else 'Not found'
   cls='exact'if status=='proven'else'number'if value else'open'
   b+='<td class="'+cls+'">'+link(f'../../index.html#n={n}&t={t}',text)+'</td>'
  b+='</tr>'
 b+='</tbody></table></div><p>Unmarked numbers are already proved maxima. “Num.” records the largest dimension attained by a checked numerical basis, including propagation to lower orders. “Not found” describes this finite search and is not a zero upper bound.</p></section>'
 b+='<section><h2>What the numerical check establishes</h2><p>For each saved basis we check orthonormality and, for every pair of basis vectors and every order up to the claimed one,</p><div class="math">'+m(r'R_t(|u_a\rangle\langle u_b|)=\delta_{ab}\frac{I_{t+1}}{t+1}.')+'</div><p>The independent implementation reads only the exported decimal coefficients and performs normalized-Dicke splitting at 80 decimal digits. It imports no optimization code. The acceptance threshold is '+m(r'10^{-10}')+' in Frobenius norm; the actual residuals are reported below. This checks the entire complex span, rather than only its individual basis vectors.</p><p>'+link('bases.json','All saved Dicke bases')+' · '+link('independent-verification.json','Independent numerical audit')+' · '+link('verify_independently.py','Reproduce the audit')+'</p></section>'
 for n in range(21,25):
  b+='<section id="n'+str(n)+'"><h2>'+str(n)+' qubits</h2><div class="scroll"><table><thead><tr><th>Order</th><th>Rigorous interval</th><th>Numerical dimension</th><th>Recorded starts</th></tr></thead><tbody>'
  for t in range(1,8):
   e=entries[n,t];interval=str(e['lower'])if e['lower']==e['upper']else str(e['lower'])+'–'+str(e['upper']);starts=e['numericalSearch']['attemptCount']
   b+='<tr><td>'+str(t)+'</td><td>'+interval+'</td><td>'+str(e['reportedDimension'] or '—')+'</td><td>'+str(starts or '—')+'</td></tr>'
  b+='</tbody></table></div><h3>Accepted numerical bases</h3><ul>'
  for c in survey['candidates']:
   if c['n']!=n:continue
   q=checks[c['id']]
   b+='<li>'+link('../../index.html#n='+str(n)+'&t='+str(c['t'])+'&subspace='+c['id'],f"Order {c['t']}, dimension {c['dimension']}")+': maximum cross-reduction error '+sci(q['maxCrossReductionFrobenius'])+'; Gram error '+sci(q['orthonormalityFrobenius'])+'.</li>'
  b+='</ul><h3>Search targets</h3><div class="scroll"><table><thead><tr><th>'+m('(t,k)')+'</th><th>Starts</th><th>Best residual</th><th>Outcome</th></tr></thead><tbody>'
  for target in survey['targets']:
   if target['n']!=n:continue
   residual=target.get('bestResidual',target.get('maxCrossReductionFrobenius'))
   b+='<tr><td>'+link(target['groupSummary'] if len(target.get('sourceFiles',[]))>1 else target['sourceFile'],f"({target['t']}, {target['dimension']})")+'</td><td>'+str(target['attemptCount'])+'</td><td>'+(sci(residual)if residual is not None else'—')+'</td><td>'+('Candidate found'if target['found']else'No accepted candidate')+'</td></tr>'
  b+='</tbody></table></div><p>Each target links to its individual record or grouped restart summary, with seeds, residual definitions and all constituent logs. The search combines unrestricted complex starts, continuation from nearby sizes, and selected structured supports. A failed target remains unresolved.</p></section>'
 b+='<section><h2>Reproducibility</h2><p>Search scripts, source snapshots and restart records are retained in the following groups. Existing exact bases and upper-bound proofs remain in the catalogue.</p><ul>'+''.join('<li>'+link(g['source'],g['group']+' search report')+'</li>'for g in survey['groups'])+'</ul><p>'+link('assemble.py','Rebuild the combined numerical export')+'</p></section>'
 css='body{margin:0;background:#f5f4ee;color:#243331;font:17px/1.65 system-ui,sans-serif}main{max-width:1120px;margin:auto;padding:36px 22px}h1,h2,h3{font-family:Georgia,serif;line-height:1.2}h1{font-size:2.5rem}a{color:#0b6356}section{border-top:1px solid #c7cfca;padding:24px 0}.scroll,.math{overflow-x:auto}table{border-collapse:collapse;width:100%;margin:12px 0}th,td{border:1px solid #ccd2cc;padding:9px 12px;text-align:center}th{background:#eaece5}.exact{background:#dfefe3}.number{background:#e4eef4}.open{background:#efede6}li{margin:.5em 0}.katex{font-size:1.08em}'
 script='document.addEventListener("DOMContentLoaded",()=>document.querySelectorAll("[data-tex]").forEach(el=>{if(window.katex)el.innerHTML=katex.renderToString(el.dataset.tex,{throwOnError:false,trust:false});}));'
 head='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Numerical AC-subspace survey · N21–24</title><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script defer src="../../vendor/katex/katex.min.js"></script><style>'+css+'</style><script>'+script+'</script></head><body><main>'
 (HERE/'index.html').write_text(head+b+'</main></body></html>\n')
if __name__=='__main__':main()
