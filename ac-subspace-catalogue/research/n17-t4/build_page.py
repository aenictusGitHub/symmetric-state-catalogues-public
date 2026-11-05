#!/usr/bin/env python3
"""Render the seventeen-qubit radical states and optimality proof offline."""
from pathlib import Path
import html
import json
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def tex(value, display=False):
    tag = 'div' if display else 'span'
    attributes = ' data-display="true"' if display else ''
    return (f'<{tag} class="math" data-tex="{html.escape(value, quote=True)}"'
            f'{attributes}></{tag}>')


def display(value):
    return tex(value, True)


def main():
    state = json.loads((HERE / 'basis.json').read_text())
    for filename in ('verification.json', 'independent-verification.json',
                     'geometry-verification.json'):
        assert json.loads((HERE / filename).read_text())['passed'], filename
    template = (ROOT / 'research/t2-spin2-inertia/index.html').read_text()
    css = re.search(r'<style>(.*?)</style>', template, re.S).group(1)
    body = '''
<p class="eyebrow">Anticoherent subspaces · Analytical states and optimality</p>
<a class="back-link" href="../../index.html#n=17&amp;t=4&amp;subspace=radical-n17-k1-t4">← Return to the catalogue</a>
<h1>Seventeen qubits at fourth order</h1>
<p class="lead">Two conjugate radical states attain fourth-order anticoherence. A saturation argument also excludes every two-dimensional fourth-order AC subspace at this size.</p>
<div class="result"><strong>Proved maximum dimension</strong>'''
    body += display(r'\boxed{K(17,4)=1}')
    body += '''<p>The explicit states below attain this maximum. Their AC order four is also globally optimal at seventeen qubits; <a href="../n17-t5/index.html">a separate exact certificate excludes every fifth-order state</a>.</p></div>
<nav class="jump-links" aria-label="On this page"><a href="#states">Analytical states</a><a href="#verification">Exact verification</a><a href="#optimality">Dimension proof</a><a href="#symmetry">Symmetry and geometry</a><a href="#sources">Sources and files</a></nav>

<section id="states"><h2>Six nonzero Dicke coefficients</h2>
<p>Use normalized Dicke states, with the superscript denoting the number of excitations. All square roots below are the positive real roots.</p>'''
    body += display(state['derivationLatex'])
    body += display(state['basis'][0]['coefficientDefinitionLatex'])
    body += '<div class="witness">' + display(state['basis'][0]['stateLatex'])
    body += display(r'|\psi_-\rangle=\overline{|\psi_+\rangle}.')
    body += '''<p>Every other Dicke coefficient is zero. Each state spans a one-dimensional AC subspace; their two-dimensional span is not AC.</p></div>
<p><a href="../../index.html#n=17&amp;t=4&amp;subspace=radical-n17-k1-t4">Inspect the exact coefficients and Majorana sphere ↗</a> · <a href="basis.json">Download the analytical coefficients ↓</a></p>
<p class="source-note">This is an independent analytical representative using the six-component support of the numerical state in <a href="../../references/ACLiege.jl">ACLiege.jl, line 180</a>, together with an additional palindromic condition. Equality or collective-rotation equivalence to that numerical state is not asserted.</p>
</section>

<section id="verification"><h2>Every reduced-state entry is checked exactly</h2>
<p>The probabilities are rational:</p>'''
    body += display(r'|A|^2=\frac{277}{1458},\qquad |B|^2=\frac{28}{243},\qquad |C|^2=\frac{142}{729}.')
    body += display(r'2\bigl(|A|^2+|B|^2+|C|^2\bigr)=1.')
    body += '''<p>For a split into four and thirteen qubits, form the amplitude matrix</p>'''
    body += display(r'(M_4)_{ab}=c_{a+b}\sqrt{\frac{\binom4a\binom{13}b}{\binom{17}{a+b}}},\qquad 0\le a\le4,\quad0\le b\le13.')
    body += display(r'\rho_4=M_4M_4^\dagger=\frac15\begin{pmatrix}1&0&0&0&0\\0&1&0&0&0\\0&0&1&0&0\\0&0&0&1&0\\0&0&0&0&1\end{pmatrix}.')
    body += '''<p>Both verifiers reconstruct all 25 entries for each conjugate branch using exact arithmetic. This proves fourth-order anticoherence, including every off-diagonal condition.</p>
<h3>Why the order is exactly four</h3><p>The fifth reduction has diagonal entries ''' + tex(r'1/6') + '''. Its only nonzero upper-triangular off-diagonal entries are</p>'''
    body += display(r'(\rho_5)_{03}=(\rho_5)_{25}=-\frac{\sqrt{8398}}{9282},\qquad (\rho_5)_{14}=\frac{\sqrt{20995}}{9282}.')
    body += display(r'\operatorname{Tr}(\rho_5^2)=\frac{5443}{32487}=\frac16+\frac{19}{21658}>\frac16.')
    body += '''<p>Thus neither representative is fifth-order AC. A <a href="../n17-t5/index.html">separate exact three-copy certificate</a> now excludes every fifth-order state at seventeen qubits, proving that these states have globally maximum AC order four.</p>
<details><summary>How the radical phases are obtained</summary><div>
<p>Palindromicity makes all odd diagonal moments vanish. Matching the zeroth, second and fourth moments fixes the three probabilities above. The support spacing eliminates shifts one, two and four. The two remaining coherences can be written as</p>'''
    body += display(r'X=48\sqrt{35}\,\overline A B=x+iy,\qquad Y=6\sqrt{10010}\,\overline B C=-D-x-2iy,')
    body += display(r'D=\frac{5680}{81},\qquad x=\frac{8(710-9q)}{243},\qquad y=\frac{8v}{243}.')
    body += display(r'2\operatorname{Re}(X+Y)+720|C|^2=0,\qquad12i\operatorname{Im}X+6i\operatorname{Im}Y=0.')
    body += '''<p>The identities follow exactly from ''' + tex(r'q^2=20995') + ' and ' + tex(r'(710-9q)^2+v^2=1628760') + '''. The complete algebra is recorded in the LaTeX derivation.</p></div></details></section>

<section id="optimality" class="section-start"><h2>Why the maximum dimension is one</h2>
<p>The argument applies to arbitrary complex AC subspaces. It makes no assumption about sparse support, palindromicity or Majorana symmetry.</p>
<div class="proof-chain"><span>Tensor packing</span><span>→ rank-sixteen deviation</span><span>→ a four-dimensional projector</span><span>→ an impossible thirty-two-qubit state</span></div>
<h3>1. Only equality in tensor packing remains</h3>
<p>For a fourth-order AC subspace of dimension ''' + tex('k') + ''', tensor operators of ranks zero, one and two generate nine orthogonal images per basis vector. The ambient symmetric space has dimension eighteen:</p>'''
    body += display(r'9k\le18\quad\Longrightarrow\quad k\le2.')
    body += '''<p>Suppose equality is attained by a rank-two projector ''' + tex('P') + '''. Let ''' + tex(r'T_{\ell m}') + ''' be Hilbert–Schmidt-orthonormal irreducible spin tensors. Fourth-order AC gives</p>'''
    body += display(r"\langle\psi_a|T_{\ell m}^{\dagger}T_{\ell' m'}|\psi_b\rangle=\frac{\delta_{ab}\delta_{\ell\ell'}\delta_{mm'}}{18}\qquad(0\le\ell,\ell'\le2).")
    body += '''<p>The eighteen tensor images therefore form a complete orthogonal basis, each with squared norm ''' + tex(r'1/18') + '''. Define the covariant completely positive map</p>'''
    body += display(r'\mathcal E_2(Z)=\sum_{\ell=0}^{2}\sum_{m=-\ell}^{\ell}T_{\ell m}ZT_{\ell m}^{\dagger}.')
    body += display(r'\mathcal E_2(P)=\frac{I_{18}}{18}.')
    body += '''<h3>2. The projector can deviate only in tensor rank sixteen</h3>
<p>On the tensor-rank-''' + tex('L') + ''' sector, the multiplier of this map is</p>'''
    body += display(r'e_{2,L}=\frac{30(x-272)(x-114)}{18\cdot323\cdot320},\qquad x=L(L+1).')
    body += '''<p>Among integer ranks ''' + tex(r'1\le L\le17') + ''', only ''' + tex(r'L=16') + ''' gives a zero. The second root, ''' + tex(r'x=114') + ''', is not of the form ''' + tex(r'L(L+1)') + '''. Consequently</p>'''
    body += display(r'P=\frac{I_{18}}9+X,\qquad X\text{ has only tensor rank }16.')
    body += '''<h3>3. Trace out one qubit</h3>
<p>The first-order AC cross conditions make the two Schmidt partners of each code basis vector orthogonal to those of the other. Therefore</p>'''
    body += display(r'Q=2\operatorname{Tr}_1P')
    body += '''<p>is a rank-four orthogonal projector in the seventeen-dimensional symmetric sixteen-qubit space. Partial trace preserves tensor rank whenever the output carries that rank, so</p>'''
    body += display(r'Q=\frac4{17}I_{17}+Y,\qquad Y\text{ is Hermitian and has only tensor rank }16.')
    body += '''<p>Its projector identity fixes the square of the nonzero deviation:</p>'''
    body += display(r'Y^2=\frac9{17}Y+\frac{52}{289}I_{17},\qquad\operatorname{Tr}(Y^2)=\frac{52}{17}.')
    body += '''<h3>4. Lift the deviation to a symmetric pure state</h3>
<p>Let ''' + tex('F') + ''' be the unitary linear part of spin-eight time reversal. Time-reversed vectorization identifies rank-sixteen operators with the maximal-spin sector of ''' + tex(r'V_8\otimes V_8') + '''. Thus</p>'''
    body += display(r'|\Omega\rangle=\frac{\operatorname{vec}(YF)}{\sqrt{\operatorname{Tr}(Y^2)}}')
    body += '''<p>is a pure symmetric thirty-two-qubit state. Its sixteen-qubit reduction is</p>'''
    body += display(r'\rho_{16}=\frac{Y^2}{\operatorname{Tr}(Y^2)}=\frac{I_{17}}{17}+\frac9{52}Y.')
    body += '''<p>Tracing one more qubit eliminates the rank-sixteen component, because the output operator space has ranks only up to fifteen. It follows that</p>'''
    body += display(r'\rho_{15}=\frac{I_{16}}{16}.')
    body += '''<p>The lifted state would therefore be fifteenth-order AC. But the tensor-Gram bound for such a pure state requires</p>'''
    body += display(r'\sum_{\ell=0}^{7}(2\ell+1)=64\le32+1=33,')
    body += '''<p>which is impossible. The hypothetical rank-two projector does not exist. Together with the explicit radical state, this proves</p>'''
    body += '<div class="certificate"><p class="eyebrow">Sharp result</p>' + display(r'\boxed{K(17,4)=1}.')
    body += '''<p>The independent verifier reconstructs all channel multipliers and checks the time-reversal and partial-trace identities used in this proof.</p></div></section>

<section id="symmetry"><h2>Two distinct rotation families</h2>
<p>Both branches have full proper rotational symmetry group ''' + tex(r'D_3') + '''. The support gives a threefold axis, and palindromicity gives a perpendicular half-turn. The support spacing excludes a higher axial order. Since the state is fourth-order but not fifth-order AC, its rank-five component is nonzero. A tetrahedral group has no rank-five invariant:</p>'''
    body += display(r'\dim\operatorname{Inv}_{T}(L=5)=\frac{11-8-3}{12}=0.')
    body += '''<p>This excludes tetrahedral symmetry and the octahedral and icosahedral groups that contain it. The constellation is not confined to an axis, so its symmetry is not continuous. The proper group is exactly ''' + tex(r'D_3') + ''', with a unique threefold axis.</p>
<p>Any rotation relating the two branches would preserve that axis up to sign. Such rotations and a global phase preserve</p>'''
    body += display(r'R=\frac{c_1c_{16}}{c_4c_{13}}=\frac{A^2}{B^2}.')
    body += '''<p>This ratio has a nonzero imaginary part. Conjugation replaces it by its distinct complex conjugate, so the two branches cannot be related by a collective rotation. These are two exact rotation families; no classification of all fourth-order states at seventeen qubits is claimed.</p>
<h3>Seventeen distinct Majorana points</h3>
<p>The exact Majorana polynomial reduces to a quintic:</p>'''
    body += display(r'P(z)=-\sqrt{17}A\,z\,Q(z^3),\qquad Q(y)=1+Ly+My^2-My^3-Ly^4-y^5.')
    body += display(r'L=-\frac{2(u+w)}{277},\qquad M=-\frac{16(2130+u+2w)(u+w)}{1628760},')
    body += display(r'u=710-9q,\qquad w=iv.')
    body += '''<p>Exact arithmetic in ''' + tex(r'\mathbb Q(w)') + ''' gives ''' + tex(r"\gcd(Q,Q')=1") + ''' and ''' + tex(r'Q(0)=1') + '''. There are one simple zero, fifteen simple nonzero roots, and one simple root at infinity. The conjugate branch has the same multiplicities. Coordinates in the viewer are numerical renderings of these certified simple roots; reconstructing the Dicke state from them has error below ''' + tex(r'2\times10^{-15}') + '''.</p></section>

<section id="sources"><h2>Source and reproducibility</h2>
<p>The six-component numerical source appears in John Martin’s <a href="../../references/ACLiege.jl">ACLiege.jl, line 180</a>. The analytical state here adds a palindromic constraint and is an independent representative. Exact reduced-matrix identities establish its existence without relying on decimal agreement with the source.</p>
<div class="downloads">
<div class="download"><a href="proof.tex">Complete LaTeX derivation ↓</a><span>Coefficients, full optimality argument, symmetry and exact multiplicities.</span></div>
<div class="download"><a href="basis.json">Analytical coefficients ↓</a><span>Exact radicals, numerical values, conventions and source distinction.</span></div>
<div class="download"><a href="verify.py">Exact state verifier ↓</a><span><a href="verification.json">Read its verification report.</a></span></div>
<div class="download"><a href="verify_independently.py">Independent verifier ↓</a><span><a href="independent-verification.json">State and upper-bound verification report.</a></span></div>
<div class="download"><a href="verify_geometry.py">Exact geometry verifier ↓</a><span><a href="geometry.json">Sphere data</a> · <a href="geometry-verification.json">geometry report</a>.</span></div>
<div class="download"><a href="build_page.py">Reproduce this page ↓</a><span>The page and its mathematical notation work entirely offline.</span></div>
</div></section>
<footer><a href="../../index.html#n=17&amp;t=4&amp;subspace=radical-n17-k1-t4">Seventeen-qubit catalogue entry ↗</a><a href="../general-bounds/index.html">General AC-subspace bounds ↗</a></footer>
'''
    # Check renderer input while preserving standalone local mathematical fonts.
    head = '''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="Exact radical fourth-order anticoherent states at seventeen qubits and a proof that the maximal AC-subspace dimension is one."><title>Seventeen qubits: exact fourth-order AC states and optimal dimension</title><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script defer src="../../vendor/katex/katex.min.js"></script><style>'''
    script = '''</style><script>document.addEventListener("DOMContentLoaded",()=>document.querySelectorAll("[data-tex]").forEach(el=>{if(window.katex)el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==="true",throwOnError:false,strict:"ignore",trust:false});}));</script></head><body><main>'''
    (HERE / 'index.html').write_text(head + css + script + body + '</main></body></html>\n')


if __name__ == '__main__':
    main()
