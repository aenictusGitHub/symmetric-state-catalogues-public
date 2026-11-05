from pathlib import Path
import html,json
HERE=Path(__file__).resolve().parent
def m(t):return '<div class="math" data-display="true" data-tex="'+html.escape(t,quote=True)+'"></div>'
def links(rows):return '<ul>'+''.join('<li><a href="'+html.escape(u,quote=True)+'">'+html.escape(label)+'</a></li>' for u,label in rows)+'</ul>'
def save(title,body):
 css='body{margin:0;background:#f5f4ee;color:#243331;font:18px/1.65 system-ui,sans-serif}main{max-width:1000px;margin:0 auto;padding:40px 24px}h1,h2{font-family:Georgia,serif;line-height:1.2}h1{font-size:2.4rem}a{color:#0b6356}section{border-top:1px solid #c7cfca;padding:25px 0}.math{overflow-x:auto;padding:12px 0}.katex{font-size:1.1em}li{margin:.4em 0}.verdict{background:#dfefe3;padding:20px;border-left:4px solid #0b6356}'
 script='document.addEventListener("DOMContentLoaded",()=>document.querySelectorAll("[data-tex]").forEach(el=>{if(window.katex)el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==="true",throwOnError:false,strict:"ignore",trust:false});}));'
 (HERE/'index.html').write_text('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script defer src="../../vendor/katex/katex.min.js"></script><style>'+css+'</style><script>'+script+'</script></head><body><main>'+body+'</main></body></html>\n')
N=19
def main():
 d=json.loads((HERE/'basis.json').read_text());a=json.loads((HERE/'independent-verification.json').read_text());assert a['passed'] and a['fullCrossReductionEntriesChecked']==64
 b='<a href="../../index.html#n='+str(N)+'&amp;t=3&amp;subspace='+d['id']+'">← AC-subspace catalogue</a><h1>Third-order maximum on '+str(N)+' qubits</h1>'
 b+='<div class="verdict">'+m('K('+str(N)+',3)=2')+'<p>An exact analytical plane attains the general upper bound. Every normalized complex superposition of its two basis vectors is third-order anticoherent.</p></div>'
 b+='<section><h2>Explicit analytical basis</h2><p>The following coefficients refer to normalized Dicke vectors. The phase definitions use positive real square roots.</p>'
 for row in d['basis']:b+=m(row['stateLatex'])
 definitions=list(dict.fromkeys(row['coefficientDefinitionLatex'] for row in d['basis']))
 b+='<p>The coefficients of both vectors are defined by the following common recurrence:</p>'
 for definition in definitions:b+=m(definition)
 b+='<p>The complete coefficient arrays are provided in machine-readable form; the construction is exact and does not depend on fitting numerical amplitudes.</p>'+links([('basis.json','Exact Dicke coefficient arrays'),('proof.tex','Complete analytical construction and proof')])+'</section>'
 b+='<section><h2>Why the whole plane is AC3</h2><p>Within each support, distinct Dicke levels differ by at least four. The diagonal probability distributions have exactly the uniform Dicke moments through degree three. The remaining off-diagonal cross terms vanish by the explicit phase closure identities.</p>'+m(r'R_3(|u_a\rangle\langle u_b|)=\delta_{ab}\frac{I_4}{4},\qquad a,b\in\{1,2\}.')+'<p>An independent verifier checks all 64 entries of the four ordered reduced matrices, including the off-diagonal operators between the basis vectors. Orthonormality follows as well.</p>'+links([('verify_independently.py','Independent exact verifier'),('independent-verification.json','Complete full-span verification report'),('verification.json','Exact physical-spin checks')])+'</section>'
 b+='<section><h2>Why dimension three is impossible</h2><p>The general auxiliary-spin-five-halves theorem applies to every third-order AC subspace at every qubit number at least nineteen:</p>'+m(r'6K(N,t)\le N-4,\qquad N\ge19,\quad t\ge3.')+m('6K('+str(N)+',3)\le '+str(N-4)+r'\quad\Longrightarrow\quad K('+str(N)+',3)\le2.')+'<p>The explicit plane meets this upper bound. The exclusion is unrestricted: no support, reality or constellation-symmetry assumption is made.</p>'+links([('../t3-spin5half-inertia/index.html','General theorem and proof'),('../t3-spin5half-inertia/independent-verification.json','Independent exact spectrum and trace checks')])+'</section>'
 b+='<section><h2>Certified Majorana representation</h2><p>The catalogue visualizes both basis vectors. Exact rational root disks certify their Majorana multiplicities; this geometry is separate from the full-subspace AC proof.</p>'+links([('geometry-verification.json','Majorana geometry certificate')])+'</section>'
 save('Third-order maximum on '+str(N)+' qubits',b)
if __name__=='__main__':main()
