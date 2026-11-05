"""Render the offline proof page from the independently checked exact certificate."""
from pathlib import Path
from html import escape
from fractions import Fraction as F
import json,hashlib
ROOT=Path(__file__).resolve().parent
checks=json.loads((ROOT/'independent-verification.json').read_text());cert=json.loads((ROOT/'sos-certificate.json').read_text());config=json.loads((ROOT/'config.json').read_text())
assert checks['fullPolynomialIdentityVerified'] and checks['attainmentVerified']
assert checks['certificateSha256']==hashlib.sha256((ROOT/'sos-certificate.json').read_bytes()).hexdigest()
k=config['k'];mode=config['mode'];word='maximum'if mode=='max'else'minimum';degree=checks['degree'];mult=(degree-4)//2;B=checks['positiveDefiniteGramBlocks'];bound=F(config['bound']);state=config['stateLatex']
def tex(t):return '<div class="proof-equation" data-tex="'+escape(t,quote=True)+'" data-display="true"></div>'
def frac(v):
 f=F(v);return str(f.numerator) if f.denominator==1 else ('-' if f<0 else '')+r'\frac{'+str(abs(f.numerator))+'}{'+str(f.denominator)+'}'
def section(id,title,s):return f'<section class="reference-section" id="{id}"><h2>{title}</h2>'+s+'</section>'
url=f'../../index.html#n=12&amp;k={k}&amp;mode={mode}'
body=f'''<section class="intro"><div class="intro-copy"><p class="eyebrow">Exact rational certificate · 4 November 2026</p><h1>Twelve qubits:<br><em>the {word} of S<sub>{k}</sub></em></h1><p class="intro-description">An exact sum-of-squares identity proves the sharp bound for every complex pure symmetric twelve-qubit state.</p></div><div class="intro-aside"><span class="section-number">12</span><p>Proved global {word}<br>Exact rational verification</p><a href="{url}">Open catalogue entry ↗</a></div></section>'''
s='<p>For every normalized pure symmetric twelve-qubit state,</p>'+tex(r'\boxed{S_'+str(k)+(r'\le'if mode=='max'else r'\ge')+frac(bound)+'.}')
s+='<p>Equality is attained by</p>'+tex(r'|\psi_\star\rangle='+state+'.')
s+='<p>Its complete sector vector is verified exactly:</p>'+tex(r'(S_1,\ldots,S_{12})=\left('+','.join(frac(v)for v in checks['attainingStateSectorVector'])+r'\right).')
s+='<p>The candidate was recorded in the catalogue before this proof. The independent certificate below establishes global optimality.</p>'
body+=section('theorem','The sharp '+word,s)
s='<p>For arbitrary, possibly unnormalized amplitudes, write</p>'+tex(r'c_r=\sqrt{\binom{12}r}\,x_r,\qquad q=\langle\psi|\psi\rangle=\sum_{r=0}^{12}\binom{12}r|x_r|^2.')
s+='<p>Define the homogeneous sector polynomial</p>'+tex('P_'+str(k)+'=q^2S_'+str(k)+r'(\psi/\sqrt q).')
lhs=(frac(bound)+'q^2-P_'+str(k))if mode=='max'else('P_'+str(k)+'-'+frac(bound)+'q^2')
s+='<p>The explicit rational data give the identity</p>'+tex(r'\boxed{q^{'+str(mult)+r'}\left('+lhs+r'\right)=\sum_{B=1}^{'+str(B)+r'}\sum_{m=0}^{2L_B}\frac{\mathbf v_{B,m}^{\dagger}G_B\mathbf v_{B,m}}{\binom{2L_B}{m}}.}')
s+=f'<p>The polynomial vectors have degree {degree//2}. Each matrix is real, symmetric and positive definite. Thus the right-hand side is nonnegative for arbitrary complex amplitudes. Dividing by the positive normalization factor gives the claimed sector bound, and the displayed analytical state proves sharpness.</p>'
s+=tex(r'q^{'+str(mult+2)+r'}>0\qquad(\psi\ne0).')
s+=f'<p>The independent verifier reconstructs the left-hand side directly from reduced density matrices. It checks all {checks["allMonomialsIncludingCancellations"]} monomial coefficients, including cancellations; {checks["nonzeroMonomials"]} are nonzero in the target. All {B} matrices pass exact leading-principal-minor tests. The largest has size {checks["maxMatrixDimension"]}. LDL decomposition gives {checks["squaredForms"]} weighted absolute squares.</p>'
s+='<p><a href="basis.json" download>Explicit polynomial vectors ↓</a> · <a href="sos-certificate.json" download>Rational Gram matrices and scales ↓</a> · <a href="independent-verification.json" download>Exact verification report ↓</a></p>'
s+='<details><summary>Polynomial conventions and reconstruction</summary><p>In <code>basis.json</code>, variable indices 0–12 denote the thirteen binomial-scaled amplitudes; indices 13–25 denote their conjugates. Each monomial is a sorted list of variable indices with a rational coefficient. Divide each channel by its corresponding <code>basisScales</code> value before applying the matrix. The positive component weights are reciprocal binomial coefficients. These finite data fully specify the proof.</p><p>The optional <a href="build_basis.py">basis generator</a> uses exact binary-form transvectants and SymPy. The independent verifier uses only Python’s standard library; no numerical solver or group-theoretic routine is required.</p></details>'
body+=section('proof',f'Proof by an exact degree-{degree} identity',s)
s='<p>Collective rotations of the displayed state attain the same bound:</p>'+tex(r'|\psi\rangle=U^{\otimes{12}}|\psi_\star\rangle,\qquad U\in SU(2).')
s+='<p>These are attaining families; completeness as a classification of all equality states is not asserted.</p>'
body+=section('family','Attaining states',s)

s='<p>The unnormalized reduced matrix is</p>'+tex(r'(\rho_t)_{ab}=\sqrt{\binom ta\binom tb}\sum_{j=0}^{12-t}\binom{12-t}{j}x_{a+j}\overline{x_{b+j}}.')
s+='<p>Reduced-purity inversion yields the sector polynomial:</p>'+tex('P_'+str(k)+'='+str(__import__('math').comb(12,k))+r'\sum_{t=0}^{'+str(k)+r'}(-1)^{'+str(k)+r'-t}\binom{'+str(k)+r'}{t}2^t\operatorname{tr}(\rho_t^2).')
s+='<p>The standalone verifier constructs this expression independently of the polynomial basis. <a href="derive_target_exact.py">Separate target derivation</a> · <a href="target-exact.json" download>Expanded target polynomial ↓</a></p>'
body+=section('derivation','Independent reduced-density-matrix derivation',s)
s='<p>From this proof directory, run:</p><pre class="precision">python3 verify_certificate.py</pre><p>The verifier uses only Python’s standard library and exact fractions. It checks the entire complex-polynomial identity, exact matrix positivity, and normalized exact-state attainment including the full sector vector.</p><p><a href="verify_certificate.py">Standalone verifier</a> · <a href="independent-verification.json">Verification report and input hashes</a></p><p>Numerical optimization was used for discovery. The final certificate and all proof checks use exact rational arithmetic; this is a computer-assisted algebraic proof with no tolerance assumptions.</p>'
body+=section('verify','Reproduce the proof',s)
html=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Twelve-qubit S{k} {word} · Exact proof</title><link rel="stylesheet" href="../../styles.css"><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script src="../../vendor/katex/katex.min.js" defer></script><style>.proof-equation{{overflow-x:auto;margin:1.4rem 0;padding:.4rem 0}}.precision{{white-space:pre-wrap;overflow-wrap:anywhere;background:var(--panel,#f0f1eb);padding:1rem}}details{{margin:1.4rem 0}}summary{{cursor:pointer}}.research-page main{{max-width:1100px}}</style></head><body class="research-page"><header class="site-header"><a class="wordmark" href="{url}">← Sector-length catalogue</a><nav aria-label="Proof navigation"><a href="#theorem">Theorem</a><a href="#proof">Proof</a><a href="#family">Attainment</a><a href="#verify">Verify</a></nav></header><main>'''+body+'''</main><footer class="site-footer"><p>Extremal symmetric states</p><a href="../index.html">All proofs and remaining gaps ↗</a></footer><script>document.addEventListener('DOMContentLoaded',()=>document.querySelectorAll('[data-tex]').forEach(el=>{el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==='true',throwOnError:false,trust:false})}));</script></body></html>'''
html=html.replace('S_10','S_{10}').replace('P_10','P_{10}').replace('S_11','S_{11}').replace('P_11','P_{11}')
(ROOT/'index.html').write_text(html)
print('Built',word,'S',k,'proof page.')
