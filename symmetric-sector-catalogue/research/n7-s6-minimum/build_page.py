"""Build the exact seven-qubit sixth-sector proof page."""
from pathlib import Path
from html import escape
from fractions import Fraction
import json,hashlib
ROOT=Path(__file__).resolve().parent
checks=json.loads((ROOT/'independent-verification.json').read_text())
cert=json.loads((ROOT/'sos-certificate.json').read_text())
assert checks['fullPolynomialIdentityVerified'] and checks['positiveDefiniteGramBlocks']==10
assert checks['certificateSha256']==hashlib.sha256((ROOT/'sos-certificate.json').read_bytes()).hexdigest()
def tex(t):return '<div class="proof-equation" data-tex="'+escape(t,quote=True)+'" data-display="true"></div>'
def frac(v):
 f=Fraction(v);return str(f.numerator) if f.denominator==1 else ('-' if f<0 else '')+r'\frac{'+str(abs(f.numerator))+'}{'+str(f.denominator)+'}'
def section(id,title,body):return f'<section class="reference-section" id="{id}"><h2>{title}</h2>'+body+'</section>'
body=r'''<section class="intro"><div class="intro-copy"><p class="eyebrow">Exact rational certificate · 3 November 2026</p><h1>Seven qubits:<br><em>the minimum of S₆</em></h1><p class="intro-description">A degree-six sum-of-squares identity proves the sharp bound for every complex pure symmetric seven-qubit state.</p></div><div class="intro-aside"><span class="section-number">7</span><p>Proved global minimum<br>Exact rational verification</p><a href="../../index.html#n=7&amp;k=6&amp;mode=min">Open catalogue entry ↗</a></div></section>'''
s=r'<p>For every normalized pure state in the symmetric subspace of seven qubits,</p>'
s+=tex(r'\boxed{S_6\ge7.}')
s+=r'<p>Equality holds exactly for collective rotations of normalized generalized GHZ states:</p>'
s+=tex(r'S_6=7\quad\Longleftrightarrow\quad |\psi\rangle=U^{\otimes7}\left(\alpha|D_7^{(0)}\rangle+\beta|D_7^{(7)}\rangle\right),\quad U\in SU(2),\quad |\alpha|^2+|\beta|^2=1.')
s+=r'<p>The coefficients may be complex. This is a complete equality classification, proved below. The coherent product state displayed in the catalogue is the case <span data-tex="\alpha=1,\ \beta=0"></span>.</p>'
body+=section('theorem','The sharp minimum',s)
s=r'<p>Write arbitrary, possibly unnormalized Dicke coefficients in binomial-scaled coordinates:</p>'
s+=tex(r'c_r=\sqrt{\binom7r}\,x_r,\qquad q=\langle\psi|\psi\rangle=\sum_{r=0}^{7}\binom7r|x_r|^2.')
s+=r'<p>Let <span data-tex="P_6=q^2S_6(\psi/\sqrt q)"></span> denote the homogeneous fourth-degree sector polynomial. The explicit rational certificate verifies</p>'
s+=tex(r'\boxed{q\bigl(P_6-7q^2\bigr)=\sum_{B=1}^{10}\sum_{m=0}^{2L_B}\frac{\mathbf v_{B,m}^{\dagger}G_B\mathbf v_{B,m}}{\binom{2L_B}{m}}.}')
s+=r'<p>Every component of <span data-tex="\mathbf v_{B,m}(x,\bar x)"></span> is a cubic polynomial with rational coefficients. Every matrix <span data-tex="G_B"></span> is real, symmetric and positive definite. Thus the right-hand side is nonnegative for arbitrary complex amplitudes. Dividing by <span data-tex="q^3&gt;0"></span> gives the universal bound. The equality family proves that the bound is sharp.</p>'
s+=f'<p>The identity was checked in exact rational arithmetic against an independent reduced-density-matrix derivation: all {checks["allMonomialsIncludingCancellations"]} coefficients agree, including cancellations. All ten matrices pass exact leading-principal-minor checks. The matrices have size at most three and give {checks["squaredForms"]} weighted absolute squares after LDL decomposition.</p>'
s+=r'<p><a href="sos-certificate.json" download>Exact Gram matrices and scales ↓</a> · <a href="basis.json" download>Explicit cubic polynomials ↓</a> · <a href="independent-verification.json" download>Independent verification report ↓</a></p>'
s+=r'<details><summary>All ten rational matrices</summary><p>These matrices apply to the scaled polynomial channels specified in the certificate. Their exact positive principal minors are supplied in the verification report.</p>'
for i,b in enumerate(cert['blocks'],1):
 matrix=r'\begin{pmatrix}'+r'\\'.join('&'.join(frac(v) for v in row) for row in b['matrix'])+r'\end{pmatrix}'
 s+=tex('G_{'+str(i)+'}='+matrix+r',\qquad L_{'+str(i)+'}='+frac(str(Fraction(b['L2'],2)))+'.')
s+=r'</details><details><summary>Polynomial and normalization conventions</summary><p>In <code>basis.json</code>, indices 0–7 denote <span data-tex="x_0,\ldots,x_7"></span>, and indices 8–15 denote their complex conjugates. A monomial is a sorted list of three variable indices with an integer coefficient. Divide each channel by its listed rational <code>basisScales</code> factor before applying the matrix. The component weights are the positive reciprocal binomial coefficients in the displayed identity. These data completely specify the proof.</p><p>The optional <a href="build_basis.py">basis generator</a> uses binary-form transvectants and the exact GHZ equality family. The verifier needs only the explicit data; it does not assume that the discovery procedure found every possible certificate.</p></details>'
body+=section('proof','Proof by an exact degree-six identity',s)
s=r'<p>Positive definiteness of every Gram matrix implies that equality forces every cubic polynomial in the certificate to vanish. The following argument classifies their common zeros.</p><h3>A rank condition from the holomorphic cubics</h3><p>The 20 holomorphic cubic components span exactly the twenty determinants of the three-column submatrices of</p>'
s+=tex(r'H=\begin{pmatrix}x_0&x_1&x_2&x_3&x_4&x_5\\x_1&x_2&x_3&x_4&x_5&x_6\\x_2&x_3&x_4&x_5&x_6&x_7\end{pmatrix}.')
s+=r'<p>The equality of these spans is verified by exact rational linear algebra. Therefore equality implies <span data-tex="\operatorname{rank}H\le2"></span>. For the binary polynomial</p>'
s+=tex(r'f(X,Y)=\sum_{r=0}^{7}\binom7r x_rX^{7-r}Y^r,')
s+=r'<p>a nonzero left-kernel vector of <span data-tex="H"></span> gives a quadratic differential equation:</p>'
s+=tex(r'(u_0\partial_X^2+u_1\partial_X\partial_Y+u_2\partial_Y^2)f=0.')
s+=r'<p>Over the complex numbers this quadratic operator factors. Distinct factors can be changed to <span data-tex="\partial_U\partial_V"></span>, so homogeneity gives a sum of two seventh powers. Repeated factors give <span data-tex="\partial_V^2"></span>, hence a tangent form. Thus the only possibilities are</p>'
s+=tex(r'f=A\ell_1^7+B\ell_2^7\qquad\text{or}\qquad f=A\ell_1^7+B\ell_1^6\ell_2.')
s+=r'<p>This argument includes all boundary cases. A single seventh power is a coherent state. A noncoherent tangent form, including one with <span data-tex="A=0"></span>, is treated separately below.</p><h3>The two coherent spinors must be orthogonal</h3><p>In the noncoherent two-power case, a collective rotation puts the first linear form along <span data-tex="X"></span>, giving</p>'
s+=tex(r'f=AX^7+B(cX+dY)^7,\qquad A\ne0,\ B\ne0,\ d\ne0.')
s+=r'<p>An unscaled highest component of the mixed cubic block with spin <span data-tex="13/2"></span> is</p>'
s+=tex(r'h=(x_0x_2-x_1^2)\overline{x_5}+(x_0x_3-x_1x_2)\overline{x_6}+(x_1x_3-x_2^2)\overline{x_7}.')
s+=r'<p>Substitution into this exact polynomial gives</p>'
s+=tex(r'h=A|B|^2c^4\bar c\,d^2\bar d^{\,5}\bigl(|c|^2+|d|^2\bigr).')
s+=r'<p>Equality requires <span data-tex="h=0"></span>. Under the stated nonzero conditions this forces <span data-tex="c=0"></span>, so the two spinors are orthogonal and the state has the generalized GHZ form.</p><h3>Noncoherent tangent states are excluded</h3><p>Rotate the repeated linear form to <span data-tex="X"></span>. The tangent form then has support only at <span data-tex="x_0,x_1"></span>. Component four, with zero-based indexing, of the same unscaled mixed cubic block becomes</p>'
s+=tex(r'-5x_1^2\overline{x_1}.')
s+=r'<p>Its vanishing forces <span data-tex="x_1=0"></span>, leaving a coherent state. Finally, every normalized generalized GHZ state has <span data-tex="S_6=7"></span>, and collective rotations preserve sector lengths. This proves both necessity and sufficiency.</p><p><a href="verify_classification.py">Independent exact classification verifier</a> · <a href="classification-independent.json">Verification record</a> · <a href="equality-classification.txt">Full argument in plain text</a> · <a href="verify_equality_classification.py">Second symbolic verifier</a></p>'
body+=section('classification','All minimizers, up to collective rotations',s)
s=r'<p>For normalized states, let <span data-tex="\rho_t"></span> be the reduced state of any <span data-tex="t"></span> qubits. Symmetry makes the subset choice irrelevant. Purity inversion and complementary purities give the exact identity</p>'
s+=tex(r'S_6=7+28\left(13\operatorname{tr}\rho_1^2+20\operatorname{tr}\rho_3^2-33\operatorname{tr}\rho_2^2\right).')
s+=r'<p>The global bound is therefore equivalent to the sharp reduced-purity inequality</p>'
s+=tex(r'\boxed{13\operatorname{tr}\rho_1^2+20\operatorname{tr}\rho_3^2\ge33\operatorname{tr}\rho_2^2.}')
s+=r'<p>Every pure symmetric seven-qubit state also has <span data-tex="S_2+S_4+S_6=63"></span>. The new bound consequently proves</p>'
s+=tex(r'\boxed{S_2+S_4\le56.}')
s+=r'<p>Both corollaries are sharp. On the GHZ family, <span data-tex="(S_0,S_2,S_4,S_6)=(1,21,35,7)"></span>. These identities were independently checked coefficient by coefficient, with arbitrary complex amplitudes. <a href="verify_even_sectors.py">Even-sector identity verifier</a> · <a href="even-sector-verification.json">Verification record</a></p>'
body+=section('corollaries','Two equivalent sharp inequalities',s)
s=r'<p>In the binomial-scaled coordinates above, the unnormalized reduced density matrix is</p>'
s+=tex(r'(\rho_t)_{ab}=\sqrt{\binom ta\binom tb}\sum_{j=0}^{7-t}\binom{7-t}{j}x_{a+j}\overline{x_{b+j}}.')
s+=tex(r'P_6=7\sum_{t=0}^{6}(-1)^{6-t}\binom6t2^t\operatorname{tr}(\rho_t^2).')
s+=r'<p>For a pure state the reductions at complementary sizes have identical nonzero spectra. The verifier checks these purity identities directly and obtains</p>'
s+=tex(r'P_6=7\left(q^2+52\operatorname{tr}\rho_1^2-132\operatorname{tr}\rho_2^2+80\operatorname{tr}\rho_3^2\right).')
s+=r'<p>The entire target <span data-tex="q(P_6-7q^2)"></span> is reconstructed from these density matrices independently of the spin-coupling construction. The expansion has 626 nonzero degree-six monomials. <a href="derive_target_exact.py">Independent target derivation</a> · <a href="target-exact.json" download>Expanded polynomial ↓</a></p>'
body+=section('derivation','Reduced-density-matrix derivation',s)
s=r'''<p>From this research directory, run the following using ordinary Python 3. All three verifiers use only the standard library and exact fractions.</p><pre class="precision">python3 verify_certificate.py
python3 verify_even_sectors.py
python3 verify_classification.py</pre><p>The optional second classification verifier uses SymPy.</p><p><a href="verify_certificate.py">Standalone global proof verifier</a> · <a href="independent-verification.json" download>Verification report and input hashes ↓</a></p><p>Numerical semidefinite optimization was used for discovery. The final certificate was reconstructed with rational coefficients; all polynomial identities and matrix positivity checks are exact. This is a computer-assisted algebraic proof covering all complex symmetric pure states.</p>'''
body+=section('verify','Reproduce the proof',s)
html=r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Seven-qubit S₆ minimum · Exact proof</title><link rel="stylesheet" href="../../styles.css"><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script src="../../vendor/katex/katex.min.js" defer></script><style>.proof-equation{overflow-x:auto;margin:1.4rem 0;padding:.4rem 0}.precision{white-space:pre-wrap;overflow-wrap:anywhere;background:var(--panel,#f0f1eb);padding:1rem}details{margin:1.4rem 0}summary{cursor:pointer}.research-page main{max-width:1100px}</style></head><body class="research-page"><header class="site-header"><a class="wordmark" href="../../index.html#n=7&amp;k=6&amp;mode=min">← Sector-length catalogue</a><nav aria-label="Proof navigation"><a href="#theorem">Theorem</a><a href="#proof">Proof</a><a href="#corollaries">Corollaries</a><a href="#verify">Verify</a></nav></header><main>'''+body+r'''</main><footer class="site-footer"><p>Extremal symmetric states</p><a href="../index.html">All proofs and remaining gaps ↗</a></footer><script>document.addEventListener('DOMContentLoaded',()=>document.querySelectorAll('[data-tex]').forEach(el=>{el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==='true',throwOnError:false,trust:false})}));</script></body></html>'''
assert not any(ord(c)<32 and c not in '\n\t' for c in html)
(ROOT/'index.html').write_text(html)
print('Built N7 S6 global proof, matrices, equality family and sharp corollaries.')
