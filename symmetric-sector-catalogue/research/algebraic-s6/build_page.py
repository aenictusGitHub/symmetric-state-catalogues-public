"""Build the static, offline algebraic S6 appendix from verified expressions."""
from pathlib import Path
from html import escape
import json
import sympy as s

ROOT = Path(__file__).resolve().parent
def read(name): return json.loads((ROOT/name).read_text())
def tex(expression, display=True):
    return f'<div class="algebraic-equation" data-tex="{escape(expression, quote=True)}" data-display="{str(display).lower()}"></div>'
def poly(expression, variable, label):
    x=s.Symbol(variable); p=s.Poly(s.sympify(expression),x)
    lines=[]
    for (power,), coefficient in p.terms():
        sign='-' if coefficient<0 else '+'
        monomial=s.latex(abs(coefficient)*x**power)
        lines.append((sign,monomial))
    chunks=[]
    for i in range(0,len(lines),2):
        chunk=''.join(sign+term for sign,term in lines[i:i+2])
        if i==0: chunk=chunk.lstrip('+')
        chunks.append((label+'&=' if i==0 else '&')+chunk)
    return tex(r'\begin{aligned}'+r'\\'.join(chunks)+r'\end{aligned}')
def table(coefficients):
    rows=''.join(f'<tr><td>{r}</td><td>{float(v):.15f}</td></tr>' for r,v in coefficients.items())
    precise='\n'.join(f'c_{r} = {v}' for r,v in coefficients.items())
    return '<table class="gap-table coefficient-table"><thead><tr><th>Dicke index r</th><th>Coefficient</th></tr></thead><tbody>'+rows+'</tbody></table><details><summary>100-digit coefficients</summary><pre class="precision">'+escape(precise)+'</pre></details>'

a10=read('n10/algebraic.json'); p10=read('n10/parametrization.json'); v10=read('n10/verified.json')
e10=read('n10/energy-polynomial.json')
a11=read('n11/exact-algebraic.json'); v11=read('n11/verification.json')
a12=read('n12/exact-algebraic.json'); v12=read('n12/verification.json')
cases=[
    dict(n=10, degree=14, value=v10['sectorLength'], coefficients=dict(zip(map(str,[0,3,6,9]),v10['coefficients'])), root=v10['root'], exact='n10/algebraic.json', verifier='n10/verify.py'),
    dict(n=11, degree=11, value=v11['sectorValue100digits'], coefficients=v11['coefficients100digits'], root=v11['ratios100digits']['B'], exact='n11/exact-algebraic.json', verifier='n11/verify_exact.py'),
    dict(n=12, degree=12, value=v12['sectorValue100digits'], coefficients={str(r):v for r,v in enumerate(v12['coefficients100digits']) if v!='0'}, root=v12['root100digits'], exact='n12/exact-algebraic.json', verifier='n12/derive_n12.py')
]
for case in cases:
    case.update(status='Exact algebraic stationary state; global optimality conjectural',globalOptimalityProven=False,exactStationarityVerified=True)
    if case['n']==10:case.update(status='Exact algebraic global minimizer',globalOptimalityProven=True,globalProof='../n10-s6-minimum/index.html')
(ROOT/'coefficients.json').write_text(json.dumps(cases,indent=2)+'\n')

body=[]
for case in cases:
    n=case['n']; label='minimizer'if n==10 else'candidate'; status='Global minimum proved by an algebraic SOS certificate'if n==10 else'Exact stationary state · global optimality conjectural'; body.append(f'<section class="reference-section" id="n{n}"><p class="eyebrow">{n} qubits · defining polynomial degree {case["degree"]}</p><h2>The {n}-qubit {label}</h2><p class="muted">{status}</p>')
    body.append(tex(r'S_6\simeq '+case['value'][:31]))
    if n==10:
        body.append('<h3 id="n10-value">An exact expression for the attained value</h3><p>The decimal above is an evaluation of the following exact algebraic number. Here the real parameters and normalization are defined by the isolated-root construction below:</p>')
        body.append(tex(r's_\star=\frac{210(1+168u+66t^2-2uz^2-480ut+96utz)}{1+120u+210t^2+10uz^2}'))
        body.append('<p>Equivalently, it is the unique real root of the integer polynomial <span data-tex="R_{14}(v)"></span> below in this rational interval:</p><pre class="precision">73.04625181143282862655926527676680 &lt; s* &lt; 73.04625181143282862655926527676682</pre>')
        body.append('<p>This value has algebraic degree 14 and Galois group <span data-tex="\\mathfrak{S}_{14}"></span>. Therefore it cannot be expressed using rational numbers, arithmetic operations and nested radicals, even allowing complex intermediate values. The exact isolated-root expression is available instead. The value is attained by the state below, and <a href="../n10-s6-minimum/index.html">global minimality is now proved</a>.</p>')
        body.append('<details class="proof-method"><summary>Minimal polynomial of the sector value</summary>')
        energy=s.Poly(s.sympify(e10['polynomial']),s.Symbol('v')).clear_denoms()[1].primitive()[1]
        body.append(poly(str(energy.as_expr()),'v','R_{14}(v)'))
        body.append('<p><a href="n10/energy-polynomial.json" download>Polynomial coefficients (JSON) ↓</a> · <a href="n10/value-verification.json" download>Exact verification certificate ↓</a> · <a href="n10/verify_value.py">Reproduce the proof</a></p></details>')
        body.append('<details class="proof-method"><summary>Why there is no radical expression</summary><p>The primitive integer polynomial stays degree 14 modulo each of the following primes, and its factors are squarefree. The exact factor degrees are:</p><table class="gap-table"><thead><tr><th>Prime</th><th>Irreducible factor degrees</th></tr></thead><tbody><tr><td>59</td><td>14</td></tr><tr><td>83</td><td>1, 2, 11</td></tr><tr><td>151</td><td>1, 13</td></tr></tbody></table><p>Irreducibility modulo 59 proves irreducibility over the rationals, hence the Galois group acts transitively on the 14 roots. The <a href="https://dummit.cos.northeastern.edu/docs/fieldthy_4_galois_theory.pdf#page=42">Dedekind–Frobenius theorem</a> gives a 13-cycle from the factorization modulo 151. Any nontrivial block system would have blocks of size 2 or 7. The 13-cycle must fix every block and then every point inside it, which is impossible. Thus the action is primitive.</p><p>The factorization modulo 83 gives a permutation with disjoint cycles of lengths 1, 2 and 11. Its eleventh power is a transposition. A primitive group containing a transposition is the full symmetric group: the pairs swapped by its conjugates form a graph whose connected components are invariant blocks; primitivity makes the graph connected, and the corresponding edge transpositions generate all permutations. Consequently the Galois group is <span data-tex="\\mathfrak{S}_{14}"></span>, which is not solvable. By <a href="https://dummit.cos.northeastern.edu/docs/fieldthy_4_galois_theory.pdf#page=47">Galois’s radical-solvability criterion</a>, no root of this irreducible polynomial has an expression by radicals over the rationals.</p></details>')
        body.append('<h3>The exact attaining state</h3>')
        body.append(tex(r'|\psi_{10}\rangle=\frac1{\sqrt Q}\left(|D_{10}^{(0)}\rangle+\sqrt{120u}|D_{10}^{(3)}\rangle+\sqrt{210}\,t|D_{10}^{(6)}\rangle+\sqrt{10u}\,z|D_{10}^{(9)}\rangle\right)'))
        body.append(tex(p10['normalizationLatex']))
        body.append('<p>Choose the isolated root '+f'<span data-tex="t\simeq {case["root"][:22]}"></span> of the degree-14 polynomial below. Set:</p>')
        body.append(tex(p10['recoveryLatex']))
        body.append('<details class="proof-method"><summary>Exact polynomial and rational reconstruction</summary>')
        for expr,label in [(a10['rootPolynomial'],'P_{14}(t)'),(a10['zNumerator'],'A(t)'),(a10['zDenominator'],'B(t)')]:
            body.append(poly(str(s.sympify(expr).subs(s.Symbol('c'),s.Symbol('t'))),'t',label))
        body.append(tex(p10['sectorLengthLatex']))
        interval=a10['rootInterval']
    elif n==11:
        body.append(tex(r'|\psi_{11}\rangle=\frac1{\sqrt Q}\left(A\sqrt T|D_{11}^{(0)}\rangle+\sqrt{165}\,B|D_{11}^{(3)}\rangle+\sqrt{462T}|D_{11}^{(6)}\rangle+\sqrt{55}|D_{11}^{(9)}\rangle\right)'))
        body.append(tex(r'Q=A^2T+165B^2+462T+55'))
        body.append(f'<p>Choose the isolated root <span data-tex="B\simeq {case["root"][:22]}"></span> of the degree-11 polynomial below. Set:</p>')
        body.append(tex(r'A=\frac{N(B)}{D(B)},\qquad T=\frac{10AB^2-10AB+706B^2+440B+39}{A^2-84AB+1848B+917}.'))
        body.append('<p>The selected branch has <span data-tex="A&lt;0"></span> and <span data-tex="T&gt;0"></span>. All square roots denote positive roots.</p><details class="proof-method"><summary>Exact polynomial and rational reconstruction</summary>')
        for key,label in [('Bpolynomial','P_{11}(B)'),('AnumeratorReduced','N(B)'),('AdenominatorReduced','D(B)')]: body.append(poly(a11[key],'B',label))
        interval=a11['Binterval']
    else:
        body.append(tex(r'\begin{aligned}|\psi_{12}\rangle={}&a\bigl(|D_{12}^{(0)}\rangle-|D_{12}^{(12)}\rangle\bigr)+b\bigl(|D_{12}^{(2)}\rangle+|D_{12}^{(10)}\rangle\bigr)\\&+c\bigl(-|D_{12}^{(4)}\rangle+|D_{12}^{(8)}\rangle\bigr)+d|D_{12}^{(6)}\rangle.\end{aligned}'))
        body.append(tex(r'2a^2+2b^2+2c^2+d^2=1'))
        body.append(f'<p>Choose the isolated root <span data-tex="x\simeq {case["root"][:22]}"></span> of the degree-12 polynomial below. Recover the coefficients by:</p>')
        body.append(tex(r'y=\frac{A(x)}{B(x)},\qquad T=\frac{2RF-PG}{2PH-RG}'))
        body.append(tex(r'c=\frac1{\sqrt{P+TR}},\quad a=\sqrt{55}\,xc,\quad b=\sqrt{14}\,y\sqrt T\,c,\quad d=\sqrt T\,c'))
        body.append('<details class="proof-method"><summary>Exact polynomial and rational reconstruction</summary>')
        for key,label in [('rootPolynomial','P_{12}(x)'),('yNumerator','A(x)'),('yDenominator','B(x)')]: body.append(poly(a12[key],'x',label))
        for key in ['P','R','F','G','H']: body.append(tex(key+'='+s.latex(s.sympify(a12[key]))))
        body.append(tex(r'S_6=\frac{4FH-G^2}{4(FR^2-GPR+HP^2)}'))
        interval=a12['rootInterval']
    body.append('<p>The defining polynomial has exactly one real root in this rational interval:</p><pre class="precision">'+escape(' < root < '.join(interval))+'</pre>')
    body.append('<p>The reconstruction denominators are nonzero and the square-root branches are positive throughout certified rational intervals. The exact stationary equations reduce to zero modulo the defining polynomial.</p></details>')
    body.append('<h3>Nonzero normalized Dicke coefficients</h3>'+table(case['coefficients']))
    body.append(f'<p><a href="{case["exact"]}" download>Exact expressions (JSON) ↓</a> · <a href="{case["verifier"]}">Verification script</a> · <a href="../../index.html#n={n}&k=6&mode=min">Catalogue record ↗</a></p></section>')

html='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Algebraic S₆ states · Extremal symmetric states</title><link rel="stylesheet" href="../../styles.css"><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script src="../../vendor/katex/katex.min.js" defer></script><style>.algebraic-equation{overflow-x:auto;margin:1.3rem 0;padding:0.3rem 0}.precision{white-space:pre-wrap;overflow-wrap:anywhere;font-size:.8rem;line-height:1.7;background:var(--panel,#f0f1eb);padding:1.2rem}.coefficient-table{max-width:570px;margin:1.5rem 0}.proof-method{margin:1.6rem 0}.research-page main{max-width:1100px}.scope-note{border-left:3px solid var(--accent,#326557);padding-left:1.2rem}summary{cursor:pointer}.research-page h3{margin-top:2rem}</style></head><body class="research-page"><header class="site-header"><a class="wordmark" href="../../index.html">← Sector-length catalogue</a><nav aria-label="Algebraic states"><a href="#n10">N = 10</a><a href="#n11">N = 11</a><a href="#n12">N = 12</a><a href="#verify">Verify</a></nav></header><main><section class="intro"><div class="intro-copy"><p class="eyebrow">Exact algebraic reconstruction · 28 October 2026</p><h1>Beyond decimal<br><em>Dicke coefficients</em></h1><p class="intro-description">Three S₆ representatives reconstructed as exact algebraic states; the ten-qubit global minimum is now proved. Isolated roots of integer polynomials determine every nonzero coefficient.</p></div><div class="intro-aside"><span class="section-number">14 · 11 · 12</span><p>Degrees of the defining polynomials<br>for N = 10, 11 and 12</p><a href="coefficients.json" download>100-digit coefficients ↓</a></div></section><section class="reference-section research-method"><h2>What has been established?</h2><p class="scope-note">These are exact algebraic stationary states matching the catalogue's numerical candidates. Their stationarity is verified in the full complex symmetric-state space. The ten-qubit state is now a proved global minimizer; the eleven- and twelve-qubit candidates remain conjectural. <a href="../n10-s6-minimum/index.html">Ten-qubit global proof ↗</a></p><p>Write the normalized state as <span data-tex="|\\psi_N\\rangle=\\sum_{r=0}^{N}c_r|D_N^{(r)}\\rangle"></span>. Rescaling the amplitudes by square roots of binomial coefficients makes the constrained stationary equations rational. Exact polynomial elimination, root isolation, and rational reconstruction give the formulas below. These polynomials were derived from the equations, not fitted to decimal data.</p><p>All omitted Dicke coefficients are zero. Collective rotations and an overall phase give equivalent representatives.</p></section>'''+''.join(body)+'''<section class="reference-section research-method" id="radicals"><h2>Why an algebraic root rather than simple radicals?</h2><p>For N = 11, the defining degree-11 polynomial has Galois group S₁₁. Consequently, this representative's coefficients cannot all be expressed by radicals over the rationals. The isolated-root form is an exact analytical description.</p><details><summary>Exact Galois-group check</summary><p>The polynomial is irreducible over the rationals, so its Galois group is transitive of prime degree 11 and therefore primitive. Modulo 71 it is squarefree with irreducible factor degrees 2 and 9. By the <a href="https://dummit.cos.northeastern.edu/docs/fieldthy_4_galois_theory.pdf#page=42">Dedekind–Frobenius theorem</a>, the corresponding permutation has cycle type (2)(9), whose ninth power is a transposition. A primitive permutation group containing a transposition is the full symmetric group. Finally, <span data-tex="B=c_3/(\\sqrt3\\,c_9)"></span>, so radical expressions for all four coefficients would give a radical expression for B. The verifier checks irreducibility and the modular factorization exactly.</p><p>Here is the elementary group step: conjugate the transposition by the group and regard the resulting swapped pairs as edges. Connected components form an invariant block system. Primitivity forces the graph to be connected, and its edge transpositions generate the full symmetric group.</p></details></section><section class="reference-section research-method" id="family"><h2>An exact family at N = 11</h2><p>There is also a continuous family with the same even sector lengths. Let <span data-tex="\\Theta|D_N^{(r)}\\rangle=(-1)^r|D_N^{(N-r)}\\rangle"></span>, extended antilinearly. For odd N, time reversal squares to −1 and the two states are orthogonal. Therefore</p>'''+tex(r'|\psi(\alpha,\beta)\rangle=\alpha|\psi_{11}\rangle+\beta\Theta|\psi_{11}\rangle,\qquad |\alpha|^2+|\beta|^2=1')+'''<p>has the same individual even-weight Pauli expectations and hence the same S₂, S₄, S₆, S₈ and S₁₀. Even-weight Pauli operators commute with time reversal, and their compression to the span of this pair is a scalar matrix. <a href="n11/analytic-notes.txt">Read the derivation</a>.</p></section><section class="reference-section research-method" id="verify"><h2>Reproduce the algebraic checks</h2><p>From the website directory, with Python 3, SymPy and mpmath:</p><pre class="precision">python3 research/algebraic-s6/n10/verify.py
python3 research/algebraic-s6/n11/verify_exact.py
python3 research/algebraic-s6/n11/verify_branches.py
python3 research/algebraic-s6/n12/derive_n12.py</pre><p>The checks include polynomial irreducibility, unique isolated roots, nonzero denominators, positive square-root branches, and exact stationarity identities. Independent high-precision evaluations reproduce the amplitudes and sector values. Numerical Hessian checks support local minimality; they do not prove global optimality.</p><p><a href="../index.html#gaps">Return to certified global bounds and remaining proof gaps ↗</a></p></section></main><footer class="site-footer"><p>Extremal symmetric states</p><a href="../../index.html">Return to catalogue ↑</a></footer><script>document.addEventListener('DOMContentLoaded',()=>{document.querySelectorAll('[data-tex]').forEach(el=>{el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==='true',throwOnError:false,trust:false})})});</script></body></html>'''
(ROOT/'index.html').write_text(html.replace('\n+python3','\npython3'))
print('Built algebraic appendix and high-precision coefficient data.')
