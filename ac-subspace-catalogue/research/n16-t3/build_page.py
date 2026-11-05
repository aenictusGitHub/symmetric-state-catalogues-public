#!/usr/bin/env python3
"""Build the offline page for the exact N=16, t=3 AC plane."""
from pathlib import Path
import html
import json
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def tex(value, display=False):
    tag = 'div' if display else 'span'
    return ('<' + tag + ' class="math" data-tex="' + html.escape(value, quote=True)
            + '"' + (' data-display="true"' if display else '') + '></' + tag + '>')


def download(path, title, description):
    assert (HERE / path).is_file(), path
    return ('<div class="download"><a href="' + html.escape(path, quote=True) + '">'
            + html.escape(title) + ' ↓</a><span>' + html.escape(description) + '</span></div>')


def main():
    state = json.loads((HERE / 'basis.json').read_text())
    report = json.loads((HERE / 'verification.json').read_text())
    independent = json.loads((HERE / 'independent-verification.json').read_text())
    geometry = json.loads((HERE / 'geometry-verification.json').read_text())
    assert report['passed'] and independent['passed'] and geometry['passed']
    assert state['id'] == 'radical-n16-k2-t3' and len(state['basis']) == 2
    css = re.search(r'<style>(.*?)</style>',
                    (ROOT / 'research/t2-spin2-inertia/index.html').read_text(), re.S).group(1)
    css += '''\n.basis-number{font-size:11px;font-weight:700;color:var(--purple);letter-spacing:.1em;text-transform:uppercase;margin-bottom:8px}.basis-block+.basis-block{margin-top:22px;padding-top:18px;border-top:1px solid var(--line)}.phase-definitions{margin-top:20px;padding-top:12px;border-top:1px solid var(--line)}.witness .math{font-size:16px}.download a{line-height:1.5}.jump-links a{white-space:nowrap}.result .math{font-size:24px}@media(max-width:600px){.witness .math{font-size:13px}.result .math{font-size:21px}}\n'''
    parts = []
    add = parts.append
    add('''<p class="eyebrow">Anticoherent subspaces · Exact radical basis</p>
<a class="back-link" href="../../index.html#n=16&amp;t=3&amp;subspace=radical-n16-k2-t3">← Return to the catalogue</a>
<h1>Sixteen qubits at third order</h1>
<p class="lead">Two analytical Dicke vectors span an optimal third-order AC subspace. Every normalized complex superposition is anticoherent, including the coherences between the basis vectors.</p>
<div class="result"><strong>Proved maximum dimension</strong>''')
    add(tex(r'\boxed{K(16,3)=2}', True))
    add('''<p>The exact two-dimensional construction attains the general dimension bound.</p></div>
<nav class="jump-links" aria-label="On this page"><a href="#basis">Analytical basis</a><a href="#proof">Full-span proof</a><a href="#maximum">Optimality</a><a href="#sources">Source and downloads</a></nav>
<section id="basis"><h2>An explicit orthonormal basis</h2><p>Use normalized Dicke states, with the superscript denoting the number of excitations. All coefficients outside the displayed supports vanish.</p><div class="witness">''')
    for i, vector in enumerate(state['basis'], 1):
        add('<div class="basis-block"><p class="basis-number">Basis vector ' + str(i) + '</p>')
        add(tex(vector['stateLatex'], True))
        add('</div>')
    add('<div class="phase-definitions"><p>The three phases are</p>')
    add(tex(state['basis'][0]['coefficientDefinitionLatex'], True))
    add('''</div></div><p><a href="../../index.html#n=16&amp;t=3&amp;subspace=radical-n16-k2-t3">Inspect the exact coefficients and Majorana spheres ↗</a> · <a href="basis.json">Download all coefficients ↓</a></p>
<p>The displayed supports are disjoint. Their squared amplitudes are</p>''')
    add(tex(r'\begin{aligned}u_1:\quad&\frac1{128}(19,20,50,20,19),\\u_2:\quad&\frac1{16}(5,3,3,5).\end{aligned}', True))
    add('''<p>Each probability list sums to one, so the two vectors are orthonormal.</p></section>
<section id="proof"><h2>Why the entire span is AC</h2><p>Write ''')
    add(tex(r'\mathcal R_t'))
    add(' for the reduction to the symmetric space of ' + tex(r't') + ''' qubits. The decisive identity is</p><div class="strong-identity">''')
    add(tex(r'\mathcal R_3(|u_a\rangle\langle u_b|)=\delta_{ab}\frac{I_4}{4},\qquad a,b\in\{1,2\}.', True))
    add('''<p>Both self reductions are maximally mixed, and both ordered cross reductions vanish.</p></div><p>For arbitrary complex coefficients, it follows that</p>''')
    add(tex(r'\begin{gathered}|\psi\rangle=\alpha|u_1\rangle+\beta|u_2\rangle,\qquad |\alpha|^2+|\beta|^2=1,\\\mathcal R_3(|\psi\rangle\langle\psi|)=\frac{I_4}{4}.\end{gathered}', True))
    add('''<h3>Checking every reduced entry</h3><p>For a coefficient vector ''')
    add(tex(r'v=\sum_{r=0}^{16}v_r|D_{16}^{(r)}\rangle'))
    add(''', form the Dicke bipartition matrix</p>''')
    add(tex(r'[M_t(v)]_{a,b}=v_{a+b}\sqrt{\frac{\binom ta\binom{16-t}b}{\binom{16}{a+b}}},\qquad\begin{array}{l}0\le a\le t,\\0\le b\le16-t.\end{array}', True))
    add(tex(r'\mathcal R_t(|v\rangle\langle w|)=M_t(v)M_t(w)^\dagger.', True))
    add('''<p>Within either basis vector, occupied Dicke indices differ by multiples of four. The off-diagonal entries of each three-qubit self reduction vanish. Inserting the two probability lists gives four equal diagonal entries, each ''')
    add(tex(r'1/4'))
    add('''.</p><p>The phase definitions cancel the remaining cross terms exactly. The primary and independent verifiers reconstruct all 64 entries of the four ordered cross reductions in rational and radical arithmetic.</p>
<details><summary>Four exact phase identities complete the proof</summary><div><p>Set ''')
    add(tex(r'F=\overline f,\quad G=\overline g'))
    add('''. The only potentially nonzero cross terms are controlled by</p>''')
    add(tex(r'\begin{aligned}E_1&=7\sqrt{19}+z\bigl(\sqrt{91}FG+15\sqrt{11}F+9\sqrt{35}G\bigr),\\E_2&=\sqrt{91}FG+\sqrt{11}F+2\sqrt{35}G,\\E_3&=2\sqrt{35}F+\sqrt{11}G+\sqrt{91},\\E_4&=7\sqrt{19}FGz+9\sqrt{35}F+15\sqrt{11}G+\sqrt{91}.\end{aligned}', True))
    add('<p>Substitution of the explicit radicals gives ' + tex(r'E_1=E_2=E_3=E_4=0') + '. Hence</p>')
    add(tex(r'M_3(u_1)M_3(u_2)^\dagger=\begin{pmatrix}0&0&E_1/896&0\\0&0&0&3zE_2/448\\3zE_3/448&0&0&0\\0&zE_4/896&0&0\end{pmatrix}=0.', True))
    add('''<p>The reverse ordered cross reduction is the adjoint. This completes the full-span verification, including arbitrary complex superpositions.</p></div></details></section>
<section id="maximum"><h2>Why dimension two is optimal</h2><p>Let ''')
    add(tex(r'\mathcal H_N=\mathrm{Sym}^N(\mathbb C^2)'))
    add(''', and consider any third-order AC subspace of dimension ''')
    add(tex(r'k'))
    add('''. Choose Hilbert–Schmidt orthonormal rectangular spin-tensor maps</p>''')
    add(tex(r'A_{\ell m}:\mathcal H_{16}\longrightarrow\mathcal H_{15},\qquad\ell=\frac12,\frac32.', True))
    add('''<p>There are six maps: two of rank one half and four of rank three halves. Every product ''')
    add(tex(r"A_{\ell m}^\dagger A_{\ell' m'}"))
    add(''' has spin-tensor ranks at most three. The whole-subspace AC condition therefore gives, for any orthonormal subspace basis ''')
    add(tex(r'v_1,\ldots,v_k'))
    add(''',</p>''')
    add(tex(r"\langle A_{\ell m}v_a,A_{\ell' m'}v_b\rangle=\frac{\delta_{ab}\delta_{\ell\ell'}\delta_{mm'}}{17}.", True))
    add('''<p>Thus six nonzero orthogonal images per basis vector must fit in the 16-dimensional output space:</p>''')
    add(tex(r'6k\le16\quad\Longrightarrow\quad k\le2.', True))
    add('''<p>The analytical plane supplies the matching lower bound. This upper bound applies to arbitrary complex subspaces, with no assumption on their supports or symmetries. <a href="../general-bounds/index.html">Read the general rectangular-tensor theorem ↗</a></p></section>
<section id="order"><h2>The two displayed states have exact AC order three</h2><p>The fourth-reduction purities are</p>''')
    add(tex(r'\begin{aligned}\operatorname{Tr}\!\left[\mathcal R_4(|u_1\rangle\langle u_1|)^2\right]&=\frac{3457}{16562}>\frac15,\\\operatorname{Tr}\!\left[\mathcal R_4(|u_2\rangle\langle u_2|)^2\right]&=\frac{3433}{16562}>\frac15.\end{aligned}', True))
    add('''<p>Neither representative is fourth-order AC, and the plane is not a fourth-order AC subspace. These two purity checks do not determine whether some other superposition in the plane has higher order.</p></section>
<section id="sources"><h2>Source and analytical reconstruction</h2><p>The starting point was John Martin’s numerical example ''')
    add(tex(r'\Pi_{16,2,3}'))
    add(''' in <a href="../../references/ACLiege.jl">ACLiege.jl, line 901</a>. Its sparse supports reveal the simple rational probability lists above. Fixing those probabilities reduces the remaining full-span conditions to phase equations; the explicit ''')
    add(tex(r'z,f,g'))
    add(''' solve them exactly.</p><p>A common collective rotation about the third axis and independent phases on the two basis vectors match the supplied numerical example:</p>''')
    add(tex(r'\begin{gathered}u_a^{\mathrm{source}}\simeq e^{i\chi_a}e^{-i\beta J_z}u_a,\\\beta\simeq-0.7306438994490132,\qquad J_z|D_{16}^{(r)}\rangle=(8-r)|D_{16}^{(r)}\rangle.\end{gathered}', True))
    add('''<p>The largest Dicke-amplitude errors are ''')
    add(tex(r'3.32\times10^{-10}'))
    add(' and ' + tex(r'5.10\times10^{-10}'))
    add('''. This decimal comparison identifies the reconstruction; the exact reduced-operator identities prove existence independently.</p>
<h3>Majorana representations</h3><p>The analytical coefficients also determine both Majorana constellations. The first vector has sixteen distinct points. The second has two double points and twelve simple points:</p>''')
    add(tex(r'\lambda(u_1)=(1^{16}),\qquad\lambda(u_2)=(2^2,1^{12}).', True))
    add('''<p>Exact polynomial square-free checks certify these multiplicities. The numerical sphere coordinates reconstruct the corresponding state rays with errors below ''')
    add(tex(r'8.4\times10^{-16}'))
    add('''.</p><h3>Proof and reproducibility</h3><div class="downloads">''')
    files = [
        ('basis.json', 'Analytical basis', 'Every exact and numerical Dicke coefficient, phase definition and source comparison.'),
        ('proof.tex', 'Complete LaTeX proof', 'Explicit cross identities and a self-contained dimension argument.'),
        ('verify.py', 'Primary exact verifier', 'Reconstructs the 64 ordered cross entries and all four phase identities.'),
        ('verification.json', 'Exact verification report', 'Saved rational and radical checks, purities and numerical source comparison.'),
        ('verify_independently.py', 'Independent verifier', 'Rebuilds the basis expressions and checks the full complex span.'),
        ('independent-verification.json', 'Independent report', 'Exact coefficient association and all ordered cross reductions.'),
        ('verify_geometry.py', 'Majorana verifier', 'Checks polynomial multiplicities and reconstruction of the state rays.'),
        ('geometry.json', 'Majorana coordinates', 'Sphere coordinates for the two analytical basis representatives.'),
        ('geometry-verification.json', 'Geometry report', 'Exact multiplicity certificates and numerical reconstruction residuals.'),
        ('build_page.py', 'Page generator', 'Rebuilds this page using the local catalogue assets.'),
    ]
    for item in files:
        add(download(*item))
    add('''</div><p class="source-note">The exact formulas reconstruct one supplied numerical plane. They prove optimal dimension without classifying every optimal plane. No claim of literature priority is made.</p></section>
<footer><span>Anticoherent subspaces · 30 October 2026</span><a href="../../index.html#n=16&amp;t=3&amp;subspace=radical-n16-k2-t3">Return to the catalogue ↗</a><a href="#top">Back to top ↑</a></footer>''')
    head = '''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="An explicit radical Dicke basis for an optimal two-dimensional third-order anticoherent subspace of sixteen symmetric qubits."><title>Sixteen qubits: an analytical third-order AC plane</title><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script defer src="../../vendor/katex/katex.min.js"></script><style>''' + css + '''</style><script>document.addEventListener("DOMContentLoaded",()=>document.querySelectorAll("[data-tex]").forEach(el=>{if(window.katex)el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==="true",throwOnError:false,strict:"ignore",trust:false});}));</script></head><body><main id="top">'''
    (HERE / 'index.html').write_text(head + '\n'.join(parts) + '</main></body></html>\n')


if __name__ == '__main__':
    main()
