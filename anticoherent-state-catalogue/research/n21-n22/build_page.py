"""Build the research page from independently verified attainment and exclusions."""
from pathlib import Path
import json,html
here=Path(__file__).resolve().parent
states=json.loads((here/'states.json').read_text())['states']
for sub in ('n21-t5','n22-t6'):
    report=json.loads((here/sub/'independent-verification.json').read_text())
    assert report['passed'] and (report.get('arbitraryComplexStatesCovered') or report.get('arbitraryComplexPureStatesCovered'))
def math(tex,display=False):
    return '<span class="math'+(' display' if display else '')+'" data-tex="'+html.escape(tex,quote=True)+'">'+html.escape(tex)+'</span>'
def links(sub,items):
    return '<div class="links">'+''.join(f'<a href="{sub}/{file}">{label} ↓</a>' for file,label in items)+'</div>'
parts=['''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>21 and 22 qubits · Proved maximum AC orders</title><link rel="stylesheet" href="../../styles.css"><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script src="../../vendor/katex/katex.min.js" defer></script><style>main{max-width:1050px;padding-block:40px}h1{font-size:clamp(32px,5vw,48px)}h2{margin-top:34px}.display{display:block;overflow-x:auto;padding:12px 0}.links{display:flex;flex-wrap:wrap;gap:16px;margin-block:20px}table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:12px;border-bottom:1px solid var(--line)}.claim{padding:18px 22px;background:var(--paper-alt,#f0efe9);border-left:3px solid var(--accent,#5e4b7c)}pre{overflow-x:auto;padding:14px;background:var(--paper-alt,#f0efe9)}</style></head><body><main><a href="../../index.html">← Pure AC-state catalogue</a><p class="eyebrow">Analytical attainment · Proven maxima</p><h1>Twenty-one and twenty-two qubits</h1><p class="claim">The maximum pure AC order is four at N=21 and five at N=22. Explicit radical states attain these orders; rigorous positivity certificates exclude the next order for every complex symmetric state.</p><table><thead><tr><th>Qubits</th><th>Maximum AC order</th><th>Universal exclusion</th></tr></thead><tbody>''']
for n,t,method in [(21,4,'Validated six-copy positivity certificate'),(22,5,'Exact rational four-copy positivity certificate')]:
    parts.append(f'<tr><th>{n}</th><td>'+math(rf't_{{\max}}({n})={t}')+f'</td><td>{method}</td></tr>')
parts.append('</tbody></table>')
for state in states:
    n,t=state['n'],state['achievedOrder'];w=state['verification']['nextOrderWitness']
    parts.append(f'<section id="n{n}"><h2>{n} qubits: maximum order {t}</h2>'+math(r'|\psi_{'+str(n)+r'}\rangle='+state['stateLatex'],True))
    parts.append('<p>All Dicke states are normalized. Exact partial traces give '+math(rf'\rho_q=I_{{q+1}}/(q+1)\quad(1\le q\le {t})')+'. The next reduction is not scalar; for example,</p>'+math(rf'(\rho_{{{t+1}}}-I_{{{t+2}}}/{t+2})_{{{w["row"]},{w["column"]}}}='+w['differenceLatex']+r'\ne0',True))
    if n==21:
        parts.append('<p>A hypothetical fifth-order state would give a positive, trace-one rotationally averaged six-copy moment whose low-multipole contractions vanish. A positive combination of these equations yields a contradiction. The complete symmetric representation and the certificate positivity are verified with outward-rounded interval arithmetic, including rigorous representation-error bounds. The verifier uses direct factorial Clebsch–Gordan contractions; floating eigenvalues are not used as proof bounds.</p><p>This closes fifth-order existence at N=21: no such pure state exists, and the displayed fourth-order state is globally optimal.</p>')
        parts.append(links('n21-t5',[('proof.tex','Complete LaTeX proof'),('certificate-bundle.zip','Complete certificate and verifier bundle'),('certificate.json','Certificate metadata'),('verify_certificate.py','Independent validated verifier'),('independent-verification.json','Validated verification report')]))
        parts.append('<p>To reproduce the proof, unpack the complete bundle, install NumPy and python-flint, and run:</p><pre><code>python3 verify_certificate.py</code></pre>')
    else:
        parts.append('<p>A hypothetical sixth-order state gives a positive, trace-one four-copy moment obeying the vanishing-multipole contraction equations. An exact rational separating identity expresses its negative trace as a sum of the vanishing equations and nonnegative matrix traces. Independent exact reconstruction gives a contradiction, proving maximum order five.</p><p>The state also spans an optimal fifth-order AC line: '+math(r'K(22,5)=1')+'. The general rectangular-tensor bound gives '+math(r'K(22,5)\le\lfloor22/12\rfloor=1')+'.</p>')
        parts.append(links('n22-t6',[('proof.tex','Complete LaTeX proof'),('certificate.json','Exact rational certificate'),('verify_independently.py','Portable independent verifier'),('independent-verification.json','Exact verification report')]))
    parts.append(f'<p><a href="../../index.html#n={n}&amp;t={t}&amp;state={state["id"]}">Open the attaining state and its Majorana sphere ↗</a></p></section>')
parts.append(links('.', [('states.json','Exact Dicke components'),('witnesses.tex','Analytical witness derivation'),('build_witnesses.py','Reproduce the constructions'),('verify_witnesses.py','Independent exact witness verifier'),('independent-verification.json','Witness verification report')]))
parts.append('<h2>Consequences and auxiliary certificates</h2><p>Every higher AC order is excluded by nesting. Since every nonzero vector of an AC subspace is an AC state, the pure-state exclusions also give '+math(r'K(21,t)=0\ (t\ge5),\qquad K(22,t)=0\ (t\ge6)')+'.</p><p>The earlier <a href="n21-t6/proof.tex">three-copy sixth-order exclusion at N=21</a> and <a href="n22-t7/proof.tex">three-copy seventh-order exclusion at N=22</a> are retained as auxiliary proofs. The stronger next-order exclusions above establish the maxima.</p></main><script>document.addEventListener("DOMContentLoaded",()=>{for(const el of document.querySelectorAll("[data-tex]"))katex.render(el.dataset.tex,el,{displayMode:el.classList.contains("display"),throwOnError:true});});</script></body></html>')
(here/'index.html').write_text(''.join(parts))
print('Built proved maximum-order page: N21=4, N22=5.')
