#!/usr/bin/env python3
"""Build the note on rigorous upper-bound consequences of nesting."""
from pathlib import Path
import html
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def tex(value, display=True):
    tag = "div" if display else "span"
    return (f'<{tag} class="math" data-display="{str(display).lower()}" '
            f'data-tex="{html.escape(value, quote=True)}"></{tag}>')


def inline(value):
    return tex(value, False)


def p(value):
    return '<p>' + value + '</p>'


def main():
    reference = (ROOT / 'research/t2-spin2-inertia/index.html').read_text()
    css = re.search(r'<style>(.*?)</style>', reference, re.S).group(1)
    css = '\n'.join(line.rstrip() for line in css.splitlines())
    css += '\n.files{display:grid;gap:.65rem}.math[data-display=false]{font-size:1em}.result .math{font-size:20px}.section-start{scroll-margin-top:24px}'
    b = '<p class="eyebrow">Anticoherent subspaces · General upper bounds</p>'
    b += '<a class="back-link" href="../t2-nested-t1/index.html">← Nesting second-order subspaces inside first-order parents</a>'
    b += '<h1>Upper bounds from nesting: what is proved</h1>'
    b += p('<strong>No stronger unrestricted upper bound is established here.</strong> Nesting gives an exact dimension constraint for a specified first-order parent. Creation and annihilation give a different, structured enlargement that recovers both existing general inequalities. The resulting catalogue bounds are unchanged.')
    b += '<div class="result"><strong>Current general upper bound</strong>'
    b += tex(r'K(N,2)\le\min\!\left\{\left\lfloor\frac{N-2}{4}\right\rfloor,\left\lfloor\frac{N+5}{5}\right\rfloor\right\},\qquad N\ge8.')
    b += p('Here ' + inline(r'K(N,2)') + ' is the largest dimension of a pure second-order anticoherent subspace of the symmetric ' + inline(r'N') + '-qubit space.') + '</div>'
    b += '<nav class="jump-links" aria-label="On this page"><a href="#same-n">First-order parents at fixed N</a><a href="#conditional">Creation and annihilation</a><a href="#interlacing">The two inequalities</a><a href="#files">Proofs and checks</a></nav>'

    b += '<section id="same-n" class="section-start"><h2>First-order parents at the same qubit number</h2>'
    b += p('Let ' + inline(r'\mathcal C') + ' be a second-order AC subspace of dimension ' + inline(r'k') + ', with projector ' + inline(r'P') + '. Write ' + inline(r'j=N/2') + ' and ' + inline(r'c=j(j+1)/3') + '. Polarization gives the whole-subspace identities')
    b += tex(r'PJ_aP=0,\qquad PJ_aJ_bP=c\delta_{ab}P.')
    b += p('The three spin images ' + inline(r'J_a\mathcal C') + ' are therefore mutually orthogonal and each has dimension ' + inline(r'k') + '. If an AC1 parent ' + inline(r'\mathcal W\supseteq\mathcal C') + ' has dimension ' + inline(r'd') + ', then first-order anticoherence also gives ' + inline(r'\langle w,J_a u\rangle=0') + ' for every ' + inline(r'w\in\mathcal W') + ' and ' + inline(r'u\in\mathcal C') + '. Hence')
    b += tex(r'\boxed{d+3k\le N+1.}')
    b += p('The exact first-order maximum ' + inline(r'K(N,1)=\lfloor N/3\rfloor+1') + ' does not supply the required lower bound on the containing parent’s dimension. A given second-order subspace need not extend to a maximum-dimensional first-order one.')
    b += '<div class="scope-note">' + p('<strong>Exact counterexample.</strong> The catalogue’s optimal ' + inline(r'N=18,\ k=4,\ t=2') + ' subspace has a largest first-order parent of dimension five, although ' + inline(r'K(18,1)=7') + '. <a href="../t2-nested-t1/index.html#example">The explicit complement matrices and proof</a> certify this obstruction.') + '</div>'
    b += p('If a maximum-dimensional parent is separately known to contain a particular code, the conditional conclusion is')
    b += tex(r'k\le\left\lfloor\frac{N-\lfloor N/3\rfloor}{3}\right\rfloor=\frac{2N}{9}+O(1).')
    b += p('For example, this gives ' + inline(r'k\le6') + ' at ' + inline(r'N=30') + ' under that extra assumption; it does not exclude arbitrary dimension-seven codes. Its asymptotic coefficient is also larger than the established ' + inline(r'1/5') + '. Using only the always-valid ' + inline(r'd\ge k') + ' gives ' + inline(r'4k\le N+1') + '.')
    b += '<details><summary>The exact extension problem</summary><div>'
    b += tex(r'R=I-P-\frac1c\sum_aJ_aPJ_a,\qquad\operatorname{rank}R=N+1-4k.')
    b += p('Every first-order parent has precisely the form')
    b += tex(r'\mathcal W=\mathcal C\oplus\mathcal A,\quad\mathcal A\subseteq\operatorname{ran}R,\quad P_{\mathcal A}J_aP_{\mathcal A}=0.')
    b += p('The complement has isotropic first and second trace moments. That is a statement about the mixed state proportional to ' + inline(r'R') + '; it does not make every vector in its support anticoherent. The three compressed operators must still be controlled. See the <a href="nesting-dimension-counts.tex">dimension-counting proof</a>.')
    b += '</div></details></section>'

    b += '<section id="conditional" class="section-start"><h2>Structured parents at adjacent qubit numbers</h2>'
    b += p('Let ' + inline(r'U') + ' be an orthonormal basis matrix of the original second-order subspace. In Schwinger notation, ' + inline(r'a,b') + ' are the two bosonic annihilation operators. Define')
    b += tex(r'E_+=\frac{(a^\dagger U,\ b^\dagger U)}{\sqrt{j+1}},\qquad E_-=\frac{(aU,\ bU)}{\sqrt j}.')
    b += p('Both are isometries with ' + inline(r'2k') + ' columns, into the symmetric spaces of ' + inline(r'N+1') + ' and ' + inline(r'N-1') + ' qubits, respectively. The AC2 identities fix their complete spin compressions:')
    b += tex(r'E_+^\dagger(J_x^\prime,J_y^\prime,J_z^\prime)E_+=\alpha_+(S_x,S_y,S_z)\otimes I_k,\qquad\alpha_+=\frac{N+3}{3},')
    b += tex(r'E_-^\dagger(J_x^\prime,J_y^\prime,J_z^\prime)E_-=\alpha_-(S_x,-S_y,S_z)\otimes I_k,\qquad\alpha_-=\frac{N-1}{3},')
    b += p('where ' + inline(r'S_a=\sigma_a/2') + '. A half-turn about the spinor ' + inline(r'y') + ' axis transforms the annihilation triple into ' + inline(r'-\alpha_-(S_x,S_y,S_z)') + '. The sign matters in the spectral argument.')
    b += p('<strong>These spaces are not first-order AC subspaces.</strong> Their spin compressions are the displayed nonzero scaled spinor matrices. They are structured creation and annihilation images, rather than same-space AC1 parents containing the original code.')
    b += '<details><summary>Why the compression coefficients are exact</summary><div>'
    b += p('For the first spinor label, the normalized diagonal entries follow from the exact Schwinger identities')
    b += tex(r'aJ_z^\prime a^\dagger=(j+J_z+1)(J_z+1/2),\quad a^\dagger J_z^\prime a=(j+J_z)(J_z-1/2).')
    b += p('Their compressions are respectively')
    b += tex(r'\frac{c+(j+1)/2}{j+1}=\frac{N+3}{6},\qquad\frac{c-j/2}{j}=\frac{N-1}{6}.')
    b += p('The raising cross entries follow from')
    b += tex(r'aJ_+^\prime b^\dagger=(j+1)^2-J_z^2,\qquad b^\dagger J_+^\prime a=j^2-J_z^2.')
    b += p('The other terms vanish by the polarized AC2 identities. This derivation applies to arbitrary complex subspaces, without any symmetry assumption. The <a href="conditional-parents.tex">full proof</a> specifies every block.')
    b += '</div></details></section>'

    b += '<section id="interlacing" class="section-start"><h2>Interlacing recovers both general bounds</h2>'
    b += '<h3>Annihilation and an auxiliary spin one</h3>'
    b += p('Couple physical spin ' + inline(r'J=(N-1)/2') + ' to an auxiliary spin one ' + inline(r'\mathbf T') + '. On the structured image, the compression of ' + inline(r'X=\mathbf J^\prime\cdot\mathbf T') + ' is equivalent to ' + inline(r'-\alpha_-\mathbf S\cdot\mathbf T\otimes I_k') + '. It has eigenvalue')
    b += tex(r'-\frac{\alpha_-}{2}=-\frac{N-1}{6}<-1\quad\text{with multiplicity }4k,\qquad N\ge8.')
    b += p('The ambient eigenvalues are ' + inline(r'-(J+1),-1,J') + ', with multiplicities ' + inline(r'N-2,N,N+2') + '. Cauchy interlacing allows at most ' + inline(r'N-2') + ' compressed eigenvalues strictly below ' + inline(r'-1') + '. Therefore')
    b += tex(r'\boxed{4k\le N-2\qquad(N\ge8).}')
    b += '<h3>Creation and an auxiliary spin three-halves</h3>'
    b += p('Now use physical spin ' + inline(r'J=(N+1)/2') + ' and auxiliary spin ' + inline(r's=3/2') + '. The compression is ' + inline(r'\alpha_+\mathbf S\cdot\mathbf T\otimes I_k') + ', with positive eigenvalue')
    b += tex(r'\frac{3\alpha_+}{4}=\frac{N+3}{4}\quad\text{with multiplicity }5k.')
    b += p('For ' + inline(r'N\ge2') + ', the second-highest ambient eigenvalue is ' + inline(r'(N-5)/4') + ', strictly below this compressed eigenvalue. The highest ambient total-spin sector has dimension ' + inline(r'N+5') + ', so interlacing gives')
    b += tex(r'\boxed{5k\le N+5.}')
    b += p('These are alternative nesting derivations of the <a href="../ancilla-inertia/index.html">established four-times-dimension bound</a> and the <a href="../t2-spin2-inertia/index.html">established five-times-dimension bound</a>. Their numerical values are unchanged.')
    b += '</section>'

    b += '<section id="files" class="section-start"><h2>Proofs and exact checks</h2>'
    b += p('The analytical proofs apply to arbitrary ' + inline(r'N') + ' in the stated ranges. The accompanying exact checks independently audit 418 Schwinger operator identities for ' + inline(r'2\le N\le20') + ', and reconstruct both conditional images of the exact ' + inline(r'N=18,\ k=4') + ' code, including their signs and normalization.')
    b += '<div class="files">'
    for path, label in [
        ('nesting-dimension-counts.tex', 'Same-qubit-number parent bound and its scope · LaTeX'),
        ('conditional-parents.tex', 'Creation/annihilation compression and interlacing proofs · LaTeX'),
        ('verify_conditional_parents.py', 'Exact operator and conditional-parent verifier'),
        ('conditional-parent-verification.json', 'Exact verification report'),
        ('../t2-nested-t1/extension-obstruction.tex', 'Exact eighteen-qubit obstruction to a maximum-dimensional parent'),
        ('../ancilla-inertia/index.html', 'Original proof of the bound 4k ≤ N − 2'),
        ('../t2-spin2-inertia/index.html', 'Original proof of the bound 5k ≤ N + 5'),
    ]:
        b += '<a href="' + path + '">' + label + ' ↗</a>'
    b += '</div>'
    b += p('These results delimit the basic dimension count and establish two exact conditional-image arguments. They do not prove that every nesting method is exhausted. A stronger universal constraint on the parent quadrupoles or the extension complement could still improve the bound.')
    b += '</section><footer><a href="../../index.html#t=2">Return to the AC-subspace catalogue ↗</a><span>No catalogue dimensions or upper bounds changed by this note.</span></footer>'

    head = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            '<title>Upper bounds from nesting: what is proved</title>'
            '<link rel="stylesheet" href="../../vendor/katex/katex.min.css">'
            '<script defer src="../../vendor/katex/katex.min.js"></script>'
            '<style>' + css + '</style>'
            '<script>document.addEventListener("DOMContentLoaded",()=>document.querySelectorAll("[data-tex]").forEach(el=>{if(window.katex)el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==="true",throwOnError:false,strict:"ignore",trust:false});}));</script>'
            '</head><body><main>')
    (HERE / 'index.html').write_text(head + b + '</main></body></html>\n')


if __name__ == '__main__':
    main()
