#!/usr/bin/env python3
"""Research note: exact construction by nesting AC orders, and its limits."""
from pathlib import Path
import html, re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

def math(s):
    return '<div class="math" data-display="true" data-tex="' + html.escape(s, quote=True) + '"></div>'

def p(s):
    return '<p>' + s + '</p>'

def main():
    css = re.search(r'<style>(.*?)</style>', (ROOT/'research/t2-spin2-inertia/index.html').read_text(), re.S).group(1)
    css = '\n'.join(line.rstrip() for line in css.split('\n'))
    css += '.math{overflow-x:auto;padding:.5rem 0;max-width:100%}.scope{padding:1rem 1.3rem;background:#f8f3e7;border-left:4px solid #8b6339}nav{display:flex;gap:1rem;flex-wrap:wrap;margin:2rem 0}.files{display:grid;gap:.7rem}table{border-collapse:collapse;width:100%}td,th{padding:.6rem;text-align:left;border-bottom:1px solid #ddd}.scroll{overflow:auto}'
    b = '<p class="eyebrow">Anticoherent subspaces · A nested construction</p><a href="../t2-general/index.html">← General construction prescriptions</a><h1>Build second-order AC subspaces inside first-order parents</h1>'
    b += p('Fixing an analytical first-order AC subspace removes all three first-moment constraints. The remaining problem involves five quadrupole matrices in a smaller space. This gives explicit analytical constructions, but restricting to maximum-dimensional first-order parents can miss optimal second-order subspaces.')
    b += '<div class="scope">' + p('<strong>An exact counterexample at 18 qubits.</strong> The catalogue’s optimal four-dimensional second-order subspace has a largest containing first-order AC subspace of dimension five. The unrestricted first-order maximum is seven. Thus even searching every seven-dimensional first-order parent would miss this particular optimal second-order subspace.') + '</div>'
    b += p('<strong>Further improvements:</strong> varying the analytical first-order parents gives certified dimensions five at N = 35 and six at N = 42. A refined radical family gives floor(N/8) − 1 for all N ≥ 40. <a href="../t2-nesting-bounds/index.html">Read the new lower bounds and constructions ↗</a>')
    b += '<nav aria-label="On this page"><a href="#reduction">Five remaining equations</a><a href="#recipe">Analytical pairing recipe</a><a href="#extension">Which parents can contain a code?</a><a href="#example">The eighteen-qubit counterexample</a><a href="#files">Proof files</a></nav>'
    b += '<section id="reduction"><h2>Compress the five quadrupoles into the parent</h2>'
    b += p('Let V be an orthonormal basis matrix for a d-dimensional first-order AC subspace. With the usual spin convention,')
    b += math(r'j=N/2,\qquad c=\frac{j(j+1)}3,\qquad V^\dagger V=I_d,\quad V^\dagger J_aV=0.')
    b += p('Write the desired k-dimensional basis as U = VX. It is second-order AC exactly when')
    b += math(r'X^\dagger X=I_k,\qquad X^\dagger A_{ab}X=0,\qquad A_{ab}=V^\dagger\!\left(\frac{J_aJ_b+J_bJ_a}{2}-c\delta_{ab}I\right)V.')
    b += p('The symmetric tensor has five independent components, since its trace is zero by the spin Casimir identity. The unknown basis has d rows instead of N + 1, and first-order anticoherence is automatic. Every equation is a whole-subspace compression identity, including cross terms between distinct basis vectors.')
    b += math(r'\{A_{ab}\}_{\rm independent}:\quad Q_{zz},\quad Q_{xx}-Q_{yy},\quad Q_{xy},\quad Q_{xz},\quad Q_{yz}.')
    b += p('For any real linear combination A of these five Hermitian matrices, its positive, negative, and zero eigenvalue counts give the necessary bound')
    b += math(r'k\le n_0(A)+\min\{n_+(A),n_-(A)\}.')
    b += p('This tests a specified parent. It does not give an unrestricted bound unless every relevant parent is covered. The catalogue’s exact maximum-dimensional first-order families are constructions, rather than a classification of all first-order subspaces.') + '</section>'
    b += '<section id="recipe"><h2>A pairing prescription with radical coefficients</h2>'
    b += p('Take one vector from each separated reflection module:')
    b += math(r'v_x=\sqrt{\frac{x+1}{2x+1}}|j,m=x\rangle+\sqrt{\frac{x}{2x+1}}|j,m=-x-1\rangle.')
    b += p('Require physical magnetic levels and separation of at least three between every two occupied levels. For example, positive x values at least one and separated by at least three meet the separation condition. All shift-one and shift-two spin terms then vanish throughout the span. Each vector has zero first moment and second moment')
    b += math(r'\langle v_x|J_z|v_x\rangle=0,\qquad h_x=\langle v_x|J_z^2|v_x\rangle=x(x+1).')
    b += p('Only the diagonal quadratic form diag(hₓ − c) remains. Pair distinct low and high modules, with no module reused:')
    b += math(r'u_{xy}=\sqrt{\frac{h_y-c}{h_y-h_x}}\,v_x+\sqrt{\frac{c-h_x}{h_y-h_x}}\,v_y,\qquad h_x<c<h_y.')
    b += p('Retain any module with hₓ = c on its own. The resulting basis is orthonormal and its entire complex span is second-order AC. It attains the exact maximum within this restricted one-vector-per-module space:')
    b += math(r'k_{\rm restricted}=\#\{h_x=c\}+\min\bigl(\#\{h_x<c\},\#\{h_x>c\}\bigr).')
    b += p('At N = 18, the modules x = 2, 5, 8 have second moments 6, 30, 72, while c = 30. Thus an exact two-dimensional example is')
    b += math(r'u_0=v_5,\qquad u_1=\sqrt{7/11}\,v_2+\sqrt{4/11}\,v_8.')
    b += p('In normalized Dicke notation Dᵣ = |D₁₈⁽ʳ⁾⟩, the coefficients are')
    b += math(r'\begin{aligned}u_0&=\sqrt{\frac6{11}}D_4+\sqrt{\frac5{11}}D_{15},\\u_1&=\frac6{\sqrt{187}}D_1+\sqrt{\frac{21}{55}}D_7+\sqrt{\frac{14}{55}}D_{12}+\frac{4\sqrt2}{\sqrt{187}}D_{18}.\end{aligned}')
    b += p('All 36 entries of the four ordered two-qubit cross reductions are verified exactly. This demonstrates the nested construction; it does not improve the known unrestricted optimum K(18,2) = 4.')
    b += '<h3>It is optimal inside the standard first-order parent</h3>'
    b += p('For the catalogue’s seven-dimensional reflection-pair parent at N = 18, compress the traceless quadrupole')
    b += math(r'A=J_z^2-30I+(J_x^2-J_y^2)+\{J_x,J_z\}.')
    b += p('Its compression has exactly five positive and two negative eigenvalues. Every second-order child has zero compression of A, so its dimension is at most two. The displayed pair attains two. At N = 24, the analogous witness')
    b += math(r'A=J_z^2-52I+2(J_x^2-J_y^2)+\{J_x,J_z\}')
    b += p('has inertia (6 positive, 3 negative) inside the nine-dimensional reflection-pair parent, giving a child dimension at most three, although the unrestricted optimum is five. Both inertia counts have exact rational interval certificates for every pivot of an LDL factorization. They apply to the specified canonical parents and their collective rotations; arbitrary phase deformations are not covered by these two finite-size certificates.')
    b += p('For the catalogue’s general reflection/central-block family, the diagonal quadrupole alone gives the asymptotic restriction')
    b += math(r'k_{\rm child}\le\frac{1-1/\sqrt3}{3}\,N+O(1)\simeq0.140883\,N.')
    b += p('It is therefore too restrictive to attain the general upper bound with leading coefficient 1/5. The <a href="reflection-parent-obstructions.tex">full parent-family proof</a> gives an explicit bound and shows failure to reach the current general upper bound for every N at least 22. This diagonal-moment argument also covers module phase changes preserving their probabilities.') + '</section>'
    b += '<section id="extension"><h2>An exact criterion for enlarging a second-order code</h2>'
    b += p('Let C already be a k-dimensional second-order AC subspace with projector P. Its polarized moment conditions imply that C and its three spin images are mutually orthogonal:')
    b += math(r'PJ_aP=0,\qquad PJ_aJ_bP=c\delta_{ab}P,\qquad \mathcal C\perp J_x\mathcal C\perp J_y\mathcal C\perp J_z\mathcal C.')
    b += p('Define their orthogonal complement projector:')
    b += math(r'R=I-P-\frac1c\sum_{a=x,y,z}J_aPJ_a,\qquad\operatorname{rank}R=N+1-4k.')
    b += p('A first-order AC parent containing C has exactly the form')
    b += math(r'\mathcal W=\mathcal C\oplus\mathcal D,\qquad\mathcal D\subseteq\operatorname{ran}R,\qquad P_{\mathcal D}J_aP_{\mathcal D}=0\quad(a=x,y,z).')
    b += p('Indeed, the first-order cross conditions force every added vector to be orthogonal to the three spin images as well as to C. Conversely those orthogonalities, together with the first-order conditions inside D, suffice. In particular,')
    b += math(r'\dim\mathcal W\le N+1-3k.')
    b += p('The exact largest parent can be found by solving a new first-order zero-compression problem for the three matrices R Jₐ R on the smaller complement.') + '</section>'
    b += '<section id="example"><h2>An optimal code whose largest first-order parent has dimension five</h2>'
    b += p('Use the <a href="../n18-t2/index.html">catalogue’s exact four-dimensional N = 18 code</a>. Here j = 9, c = 30 and R has rank three. Let eᵣ denote the normalized Dicke basis. The three columns')
    b += math(r'q_+=\frac{Re_1}{\sqrt{R_{11}}},\qquad q_0=\frac{Re_0}{\sqrt{R_{00}}},\qquad q_-=\frac{Re_2}{\sqrt{R_{22}}},\qquad Q=(q_+,q_0,q_-)')
    b += p('form an orthonormal basis of the complement. Exact radical arithmetic gives')
    b += math(r'Q^\dagger J_zQ=\frac92\begin{pmatrix}1&0&0\\0&0&0\\0&0&-1\end{pmatrix},\qquad Q^\dagger J_+Q=\frac9{\sqrt2}\begin{pmatrix}0&0&0\\1&0&0\\0&-i&0\end{pmatrix}.')
    b += p('Any two-dimensional subspace with zero compression of the first matrix must contain the middle coordinate vector. An orthogonal second vector has the form (a, 0, b), with equal nonzero magnitudes of a and b. Its cross matrix element with the middle vector under the second matrix is 9a/√2, which cannot vanish. Hence the complement allows only one added first-order dimension. The line spanned by q₀ attains it.')
    b += math(r'\boxed{\max_{\mathcal W\supseteq\mathcal C,\;\mathcal W\ {\rm AC1}}\dim\mathcal W=5\;<\;K(18,1)=7.}')
    b += p('Two separate exact verifiers reconstruct the projector and these matrices. One also reconstructs the full 19 × 16 orthogonal error frame and explicitly checks the five-dimensional parent. This is a statement about this particular optimal second-order code; it does not classify all optimal codes at N = 18.')
    b += p('Nesting remains a useful construction method. A complete search must allow parents beyond the known maximum-dimensional families, including first-order parents that are maximal under inclusion at a smaller dimension.') + '</section>'
    b += '<section id="files"><h2>Exact proofs and verification</h2><div class="files">'
    for path, label in [('nested-examples.tex','Pairing construction: LaTeX proof'),('nested_examples.py','Exact nested-basis verifier'),('nested-examples.json','Explicit basis and verification'),('extension-obstruction.tex','First-order extension theorem and counterexample'),('extension_obstruction.py','Exact extension verifier'),('extension-obstruction.json','Extension verification report'),('extension_independent_audit.py','Independent error-frame verifier'),('extension-independent-audit.json','Independent audit report'),('reflection-parent-obstructions.tex','Bounds inside the analytical parent family'),('verify_reflection_inertia.py','Specified-parent inertia verifier'),('reflection-inertia-verification.json','Certified inertia counts')]:
        b += '<a href="'+path+'">'+label+' ↗</a>'
    b += '</div>' + p('Catalogue dimensions and intervals are unchanged. The construction and extension results concern the stated subspaces and parent families. The definition and polarization conventions agree with <a href="https://doi.org/10.1103/PhysRevA.111.022435">Serrano-Ensástiga, Chryssomalakos and Martin, Physical Review A 111, 022435 (2025)</a>.') + '</section><footer><a href="../../index.html#n=18&t=2&subspace=radical-n18-k4-t2">Return to the N = 18 catalogue entry ↗</a></footer>'
    head = '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Nested constructions of second-order anticoherent subspaces</title><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script defer src="../../vendor/katex/katex.min.js"></script><style>'+css+'</style><script>document.addEventListener("DOMContentLoaded",()=>document.querySelectorAll("[data-tex]").forEach(el=>{if(window.katex)el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:true,throwOnError:false,strict:"ignore",trust:false});}));</script></head><body><main>'
    (HERE/'index.html').write_text(head+b+'</main></body></html>\n')

if __name__ == '__main__':
    main()
