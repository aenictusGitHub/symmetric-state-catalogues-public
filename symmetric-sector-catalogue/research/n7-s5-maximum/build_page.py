"""Build the offline proof page for the sharp seven-qubit fifth-sector maximum."""
from pathlib import Path
from html import escape
from fractions import Fraction as F
import json,hashlib
ROOT=Path(__file__).resolve().parent
checks=json.loads((ROOT/'independent-verification.json').read_text());cert=json.loads((ROOT/'sos-certificate.json').read_text())
assert checks['fullPolynomialIdentityVerified'] and checks['positiveDefiniteGramBlocks']==15
assert checks['certificateSha256']==hashlib.sha256((ROOT/'sos-certificate.json').read_bytes()).hexdigest()
def tex(t):return '<div class="proof-equation" data-tex="'+escape(t,quote=True)+'" data-display="true"></div>'
def frac(v):
 f=F(v);return str(f.numerator) if f.denominator==1 else ('-' if f<0 else '')+r'\frac{'+str(abs(f.numerator))+'}{'+str(f.denominator)+'}'
def section(id,title,s):return f'<section class="reference-section" id="{id}"><h2>{title}</h2>'+s+'</section>'
body='''<section class="intro"><div class="intro-copy"><p class="eyebrow">Exact rational certificate · 3 November 2026</p><h1>Seven qubits:<br><em>the maximum of S₅</em></h1><p class="intro-description">A degree-six sum-of-squares identity proves the sharp bound for every complex pure symmetric seven-qubit state.</p></div><div class="intro-aside"><span class="section-number">7</span><p>Proved global maximum<br>Exact rational verification</p><a href="../../index.html#n=7&amp;k=5&amp;mode=max">Open catalogue entry ↗</a></div></section>'''
s='<p>For every normalized pure symmetric seven-qubit state,</p>'+tex(r'\boxed{S_5\le\frac{267}{7}.}')
s+='<p>The seven-qubit W state attains equality:</p>'+tex(r'|\psi_\star\rangle=|D_7^{(1)}\rangle=\frac{|1000000\rangle+|0100000\rangle+\cdots+|0000001\rangle}{\sqrt7}.')
s+='<p>Its full sector vector is verified exactly:</p>'+tex(r'(S_1,\ldots,S_7)=\left(\frac{25}{7},\frac{51}{7},\frac{125}{7},35,\frac{267}{7},\frac{145}{7},\frac{31}{7}\right).')
s+='<p>The working draft reports this attaining candidate in Table III. The independent global certificate below establishes its optimality.</p>'
body+=section('theorem','The sharp maximum',s)
s='<p>Write arbitrary, possibly unnormalized Dicke coefficients as</p>'+tex(r'c_r=\sqrt{\binom7r}\,x_r,\qquad q=\langle\psi|\psi\rangle=\sum_{r=0}^7\binom7r|x_r|^2.')
s+=r'<p>Let <span data-tex="P_5=q^2S_5(\psi/\sqrt q)"></span> be the homogeneous fourth-degree sector polynomial. The explicit rational data verify the identity</p>'+tex(r'\boxed{q\left(\frac{267}{7}q^2-P_5\right)=\sum_{B=1}^{15}\sum_{m=0}^{2L_B}\frac{\mathbf v_{B,m}^{\dagger}G_B\mathbf v_{B,m}}{\binom{2L_B}{m}}.}')
s+=r'<p>Each component of <span data-tex="\mathbf v_{B,m}(x,\bar x)"></span> is an explicit cubic polynomial with rational coefficients. All fifteen matrices <span data-tex="G_B"></span> are real, symmetric and positive definite. Every summand is therefore nonnegative. For a nonzero state, division by <span data-tex="q^3&gt;0"></span> proves the bound. The displayed state proves sharpness.</p>'
s+=f'<p>The independent verifier reconstructs the left-hand side from reduced density matrices, expands both sides with formal complex variables, and checks all {checks["allMonomialsIncludingCancellations"]} coefficients, including cancellations ({checks["nonzeroDegreeSixMonomials"]} are nonzero in the target). Matrix positivity is checked using exact leading principal minors. The matrices have size at most three and yield {checks["squaredForms"]} weighted absolute squares after LDL decomposition.</p>'
s+='<p><a href="basis.json" download>Explicit cubic polynomials ↓</a> · <a href="sos-certificate.json" download>Rational Gram matrices and scales ↓</a> · <a href="independent-verification.json" download>Exact verification record ↓</a></p><details><summary>All fifteen rational matrices</summary>'
for i,b in enumerate(cert['blocks'],1):
 mat=r'\begin{pmatrix}'+r'\\'.join('&'.join(frac(v) for v in row) for row in b['matrix'])+r'\end{pmatrix}'
 s+=tex('G_{'+str(i)+'}='+mat+r',\qquad L_{'+str(i)+'}='+frac(F(b['L2'],2))+'.')
s+='</details><details><summary>Polynomial data conventions</summary><p>In <code>basis.json</code>, variable indices 0–7 stand for the eight binomial-scaled amplitudes, and indices 8–15 for their conjugates. Each monomial is a sorted list of three indices and a rational coefficient. Divide each channel by its corresponding <code>basisScales</code> value in the certificate before applying its matrix. Component weights are reciprocal binomial coefficients. These explicit finite data completely specify the identity; the proof does not rely on the numerical discovery procedure.</p><p>The optional <a href="build_basis.py">basis generator</a> constructs binary-form transvectants and restricts them to the W-state representative. It uses SymPy. The standalone verifier needs only Python’s standard library.</p></details>'
body+=section('proof','Exact sum-of-squares proof',s)
s='<p>Collective rotations of the W state attain the same maximum:</p>'+tex(r'|\psi\rangle=U^{\otimes7}|D_7^{(1)}\rangle,\qquad U\in SU(2).')
s+='<p>This includes the six-excitation Dicke state. This proof establishes an attaining orbit; it does not claim a complete classification of all maximizers.</p>'
body+=section('family','Attaining states',s)
s='<p>Let the state be normalized, and let the one- and two-qubit reduced purities be</p>'+tex(r'R_1=\operatorname{tr}(\rho_1^2),\qquad R_2=\operatorname{tr}(\rho_2^2).')
s+='<p>Purity inversion and equality of complementary reduced purities give</p>'+tex(r'S_5=21(10R_1-8R_2-1).')
s+='<p>The proved sector bound is therefore equivalent to the sharp inequality</p>'+tex(r'\boxed{5R_1-4R_2\le\frac{69}{49}.}')
s+='<p>For the W state, the purities are</p>'+tex(r'R_1=\frac{37}{49},\qquad R_2=\frac{29}{49},')
s+='<p>which saturate the inequality. The verifier checks both the homogeneous sector identity and these exact purities.</p>'
body+=section('purities','An equivalent sharp purity inequality',s)
s='<p>In the above coordinates, the unnormalized reduced matrix is</p>'+tex(r'(\rho_t)_{ab}=\sqrt{\binom ta\binom tb}\sum_{j=0}^{7-t}\binom{7-t}{j}x_{a+j}\overline{x_{b+j}}.')
s+='<p>Reduced-purity inversion gives the sector polynomial directly:</p>'+tex(r'P_5=21\sum_{t=0}^5(-1)^{5-t}\binom5t2^t\operatorname{tr}(\rho_t^2).')
s+='<p>The verifier uses this formula independently of the transvectant basis. <a href="derive_target_exact.py">Separate target derivation</a> · <a href="target-exact.json" download>Expanded target polynomial ↓</a></p>'
body+=section('derivation','Independent density-matrix derivation',s)
s='<p>Download this proof directory and run:</p><pre class="precision">python3 verify_certificate.py</pre><p>The verifier uses only Python’s standard library and exact fractions. It checks the full complex polynomial identity, positive definiteness of every Gram matrix, the normalized W state and its complete sector vector, and the equivalent sharp reduced-purity inequality.</p><p><a href="verify_certificate.py">Standalone verifier</a> · <a href="independent-verification.json">Results and input hashes</a></p><p>Numerical optimization was used to discover the certificate. The final proof is entirely rational and has no floating-point tolerance assumptions. This is a computer-assisted algebraic proof.</p>'
body+=section('verify','Reproduce every exact check',s)
html='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Seven-qubit S₅ maximum · Exact proof</title><link rel="stylesheet" href="../../styles.css"><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script src="../../vendor/katex/katex.min.js" defer></script><style>.proof-equation{overflow-x:auto;margin:1.4rem 0;padding:.4rem 0}.precision{white-space:pre-wrap;overflow-wrap:anywhere;background:var(--panel,#f0f1eb);padding:1rem}details{margin:1.4rem 0}summary{cursor:pointer}.research-page main{max-width:1100px}</style></head><body class="research-page"><header class="site-header"><a class="wordmark" href="../../index.html#n=7&amp;k=5&amp;mode=max">← Sector-length catalogue</a><nav aria-label="Proof navigation"><a href="#theorem">Theorem</a><a href="#proof">Proof</a><a href="#family">Maximizers</a><a href="#verify">Verify</a></nav></header><main>'''+body+'''</main><footer class="site-footer"><p>Extremal symmetric states</p><a href="../index.html">All proofs and remaining gaps ↗</a></footer><script>document.addEventListener('DOMContentLoaded',()=>document.querySelectorAll('[data-tex]').forEach(el=>{el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==='true',throwOnError:false,trust:false})}));</script></body></html>'''
(ROOT/'index.html').write_text(html)
print('Built exact S5 maximum proof page.')
