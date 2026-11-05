#!/usr/bin/env python3
"""Render rigorous improvements from nested first-order AC constructions."""
from pathlib import Path
import html,json,re
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def tex(x):return '<div class="math" data-display="true" data-tex="'+html.escape(x,quote=True)+'"></div>'
def p(x):return '<p>'+x+'</p>'
def main():
    audits=[json.loads((HERE/f'n{n}-k{k}'/'independent-verification.json').read_text())for n,k in [(35,5),(42,6)]]
    assert all(x['passed']for x in audits)
    assert json.loads((HERE/'lower-fourband-strengthened-certificate.json').read_text())['passed']
    css=re.search(r'<style>(.*?)</style>',(ROOT/'research/t2-nested-t1/index.html').read_text(),re.S).group(1)
    css='\n'.join(s.rstrip()for s in css.split('\n'))
    css+=' .result{margin:1.7rem 0}.proof-links{display:flex;flex-wrap:wrap;gap:.8rem;margin:1rem 0}details{margin:1.5rem 0}summary{cursor:pointer;font-weight:600}.record{border:1px solid #deddd3;background:#fffef9;padding:1.2rem;margin:1rem 0}.scroll{overflow-x:auto}td,th{vertical-align:top}'
    body='<p class="eyebrow">Anticoherent subspaces · Improved lower bounds from nesting</p><a href="../t2-nested-t1/index.html">← First-order parents and nesting</a><h1>Nesting gives stronger second-order constructions</h1>'
    body+=p('Allowing broader first-order parents produces exact five- and six-dimensional examples at 35 and 42 qubits. A revised radical construction also improves the general four-band lower bound by two dimensions for every N at least 40.')
    body+='<div class="result"><strong>A stronger general lower bound</strong>'+tex(r'\boxed{K(N,2)\ge\left\lfloor\frac N8\right\rfloor-1\qquad(N\ge40).}')+'</div>'
    body+='<div class="scope">'+p('<strong>Scope.</strong> These are rigorous existence results. The upper bounds are unchanged, and no infinite family attaining the general upper bound is established. The finite examples have exact isolated algebraic definitions; short radical expressions for those two bases are not claimed.')+'</div>'
    body+='<nav aria-label="On this page"><a href="#improvements">Improved entries</a><a href="#parents">Broader first-order parents</a><a href="#finite">Exact algebraic bases</a><a href="#family">General radical family</a><a href="#files">Proof files</a></nav>'
    body+='<section id="improvements"><h2>Two catalogue intervals improve</h2><div class="scroll"><table><thead><tr><th>Qubits</th><th>Previous bounds</th><th>Improved bounds</th><th>Exact witness</th></tr></thead><tbody>'
    for n,k,old,upper in [(35,5,4,8),(42,6,5,9)]:
        body+=f'<tr><td>{n}</td><td>{old}–{upper}</td><td><strong>{k}–{upper}</strong></td><td><a href="../../index.html#n={n}&amp;t=2&amp;subspace=algebraic-parity-n{n}-k{k}-t2">{k}-dimensional parity construction ↗</a></td></tr>'
    body+='</tbody></table></div>'+p('Both lower witnesses are certified for every normalized complex superposition. Their basis vectors can be inspected in the catalogue, including their Majorana constellations. Existing representatives remain available.')+'</section>'
    body+='<section id="parents"><h2>A larger analytical family of first-order parents</h2>'
    body+=p('Restrict the Dicke support to even excitation number. Then every operator shifting that number by one has zero compression. In particular, the x and y first-spin conditions hold automatically. Order the remaining levels by the sign of their magnetic number:')
    body+=tex(r'J_z\big|_{\mathcal H_{\rm even}}=\begin{pmatrix}A&0\\0&-B\end{pmatrix},\qquad A,B>0.')
    body+=p('For the two cases here the positive and negative blocks have equal dimension d: d = 9 at N = 35 and d = 11 at N = 42. Given any d × d unitary matrix T, define')
    body+=tex(r'Z(T)=\begin{pmatrix}A^{-1/2}\\B^{-1/2}T\end{pmatrix},\qquad H(T)=A^{-1}+T^\dagger B^{-1}T,\qquad V(T)=Z(T)H(T)^{-1/2}.')
    body+=tex(r'V(T)^\dagger V(T)=I_d,\qquad V(T)^\dagger J_zV(T)=0,\qquad V(T)^\dagger J_xV(T)=V(T)^\dagger J_yV(T)=0.')
    body+=p('Thus every T gives an analytical first-order AC parent. Every neutral child in these balanced parity sectors extends to a parent of this form: its positive and negative components define a partial isometry in the weighted coordinates, which extends to a unitary T. This parametrizes the relevant maximum neutral planes inside the chosen parity sector, without asserting a classification of unrestricted first-order AC subspaces.')
    body+=p('Inside such a parent, the five general quadrupole conditions reduce to three nonautomatic Hermitian conditions. Equivalently, for the desired basis U, it remains to impose')
    body+=tex(r'U^\dagger J_z^2U=cI_k,\qquad U^\dagger J_+^2U=0,\qquad c=\frac{j(j+1)}3,\quad j=N/2.')
    body+=p('The mixed tensors involving one ladder step vanish by parity. Allowing the unitary T to vary avoids the restrictions of the fixed reflection-pair parents examined in the earlier note.')+'</section>'
    body+='<section id="finite"><h2>Exact algebraic bases at 35 and 42 qubits</h2>'
    body+=p('The numerical candidates were converted to square systems of polynomial equations with rational coefficients. Unselected coordinates are fixed to explicit rational numbers. The remaining coordinates are specified as the unique root in a rational box around the recorded rational centre:')
    body+=tex(r'F(x_*)=0,\qquad\|x_*-x_0\|_\infty\le10^{-50}.')
    body+=p('For a recorded rational approximate inverse R of the Jacobian, the map x ↦ x − RF(x) is proved to be a strict contraction carrying that box into itself. Every inequality in the certificate uses integer products and exact rational comparisons. The isolating box, polynomial system, fixed coordinates, and coefficient formulas together specify each basis exactly.')
    body+=tex(r'|u_a\rangle=\sum_{s=0}^{M-1}\frac{\sqrt{\binom N{2s}}}{2^{e_s}}\bigl(z_{sk+a}+i z_{Mk+sk+a}\bigr)|D_N^{(2s)}\rangle,\qquad a=0,\ldots,k-1.')
    body+=p('Here M is the number of even Dicke levels: 18 and 22 respectively. The exponents eₛ and the full real coordinate vector z are specified in the certificate; z combines the isolated root coordinates with the fixed rational coordinates. Every odd Dicke coefficient is exactly zero.')
    for n,k,d in [(35,5,9),(42,6,11)]:
        a=next(x for x in audits if x['n']==n)
        body+=f'<div class="record"><h3>{n} qubits · dimension {k}</h3>'
        body+=p(f'The square system has {a["realVariables"]} real unknowns and {a["realVariables"]} quadratic equations, with {a["fixedRationalCoordinates"]} further coordinates fixed rationally. An independent verifier reconstructs all {a["fullCrossReductionEntriesChecked"]} entries of the ordered two-qubit cross reductions from the normalized Dicke splitting formulas.')
        body+=tex(r'M_2(u_a)M_2(u_b)^\dagger=\frac{\delta_{ab}}3 I_3.')
        body+='<div class="proof-links">'
        for file,label in [('basis.json','Exact basis'),('system.json','Polynomial system'),('certificate.json','Rational root certificate'),('verification.json','Primary verification'),('independent-verification.json','Independent verification'),('geometry.json','Majorana coordinates'),('geometry-verification.json','Geometry checks'),('verify_geometry.py','Geometry verifier')]:body+=f'<a href="n{n}-k{k}/{file}">{label} ↗</a>'
        body+='</div></div>'
    body+=p('The independent cross-reduction checks certify whole subspaces, including off-diagonal pairs of basis vectors. Rounded displayed coefficients are only a rendering of the exact algebraic definitions. Neither certificate proves optimality of dimension five or six.')+'</section>'
    body+='<section id="family"><h2>A stronger family with radical coefficients</h2>'
    body+=p('Write N = 8q + ρ, with q at least five and 0 ≤ ρ ≤ 7. The improved construction has k = q − 1 rows. Choose')
    body+=tex(r'h=\begin{cases}2q,&\rho\le5,\\2q+2,&\rho\ge6,\end{cases}\qquad o=\begin{cases}2+\lfloor\rho/2\rfloor,&\rho\le5,\\\lfloor(\rho-2)/2\rfloor,&\rho\ge6.\end{cases}')
    body+=tex(r'r_{i\ell}=o+2i+\ell h,\qquad x_i=\frac{N/2-o-2i}{h},\qquad v=\frac{N(N+2)}{12h^2},\qquad0\le i<k,\quad0\le\ell\le3.')
    body+=p('The row weights are polynomial functions of x, with a small correction in remainder class five:')
    body+=tex(r'\begin{aligned}w_0&=v/4+x^2/4-35x/36+5/6,&w_1&=-v/4-x^2/4+5x/12+1/2,\\w_2&=-v/4-x^2/4+13x/12-1/2,&w_3&=v/4+x^2/4-19x/36+1/6,\end{aligned}')
    body+=tex(r'\widetilde{\mathbf w}(x)=\mathbf w(x)+\mathbf1_{\rho=5}\frac{x-3/2}{2q}\left(-\frac13,1,-1,\frac13\right).')
    body+=p('They are strictly positive and have exactly the required normalization, zero first moment, and isotropic second moment. The analytical basis has the form')
    body+=tex(r'|u_i\rangle=\sum_{\ell=0}^3\sqrt{\widetilde w_\ell(x_i)}\,\eta_{i\ell}|D_N^{(r_{i\ell})}\rangle,\qquad|\eta_{i\ell}|=1.')
    body+=p('The phases are assigned successively by closing four lengths into a quadrilateral. The complete radical recurrence is given in the proof. Nonnegative integer-polynomial coefficient certificates prove positivity and strict quadrilateral feasibility in all eight remainder classes and for every q at least five. This is an all-N proof; the supplementary numerical checks validate the implementation only.')
    body+=p('Every four-level row has two positive and two negative magnetic levels. Its mean-zero line extends analytically to a two-dimensional neutral plane. Their disjoint, same-parity direct sum is a 2k-dimensional first-order AC parent containing the constructed k-dimensional second-order subspace.')
    body+=tex(r'\boxed{K(N,2)\ge\max\!\left\{\left\lfloor\frac{N+3}{9}\right\rfloor,\left\lfloor\frac N8\right\rfloor-1\right\}\qquad(N\ge40).}')
    body+=p('The previous four-band guarantee was floor(N/8) − 3. The stronger constant raises the combined general lower bound, for example, from eleven to twelve at N = 104. Its asymptotic coefficient remains 1/8, below the current upper-bound coefficient 1/5.')
    body+='<div class="proof-links"><a href="lower-fourband-strengthened-proof.tex">Complete radical prescription and all-N proof ↗</a><a href="lower_fourband_cert_all.py">Exact polynomial verifier ↗</a><a href="lower-fourband-strengthened-certificate.json">All-N certificate ↗</a></div></section>'
    body+='<section id="files"><h2>Verification and reproducibility</h2><div class="files">'
    files=[('proof.tex','Finite algebraic constructions and proof'),('verify_parity_certificate.py','Primary finite-basis verifier'),('verify_parity_independently.py','Independent finite-basis verifier'),('certify_parity.py','Finite certificate generator'),('lower_fourband_cert_all.py','Exact all-N coefficient verifier'),('lower-fourband-strengthened-certificate.json','General construction certificate'),('lower_fourband_independent_audit.py','Independent general-construction audit'),('lower-fourband-independent-audit.json','Independent general-construction report'),('lower_fourband_render.py','Radical construction and parent rendering'),('lower-fourband-render-verification.json','Supplementary spin-compression checks'),('../t2-nested-t1/index.html','Nesting conditions and fixed-parent limitations'),('../t2-general/index.html','Earlier general constructions')]
    for file,label in files:body+=f'<a href="{file}">{label} ↗</a>'
    body+='</div>'+p('No upper-bound interval endpoint is changed by these existence constructions. Unsuccessful searches for still larger dimensions are not exclusion proofs.')+'</section><footer><a href="../../index.html">Return to the catalogue ↗</a></footer>'
    head='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Improved AC-subspace bounds from nesting</title><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script defer src="../../vendor/katex/katex.min.js"></script><style>'+css+'</style><script>document.addEventListener("DOMContentLoaded",()=>document.querySelectorAll("[data-tex]").forEach(el=>{if(window.katex)el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:true,throwOnError:false,strict:"ignore",trust:false});}));</script></head><body><main>'
    (HERE/'index.html').write_text(head+body+'</main></body></html>\n')
if __name__=='__main__':main()
