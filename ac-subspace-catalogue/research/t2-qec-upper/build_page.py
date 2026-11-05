#!/usr/bin/env python3
"""Render QEC and joint-compression upper-bound research with explicit scope."""
from pathlib import Path
import html,json,re
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def math(x,display=True):
 tag='div'if display else'span'
 return f'<{tag} class="math" data-display="{str(display).lower()}" data-tex="{html.escape(x,quote=True)}"></{tag}>'
def p(x):return '<p>'+x+'</p>'
def link(path,label):return f'<a href="{path}">{label} ↗</a>'
def main():
 audit=json.loads((HERE/'lp-feasibility-independent-audit.json').read_text());assert audit['passed']
 search=json.loads((HERE/'matrix-inertia-search-report.json').read_text())
 css=re.search(r'<style>(.*?)</style>',(ROOT/'research/t2-nesting-upper/index.html').read_text(),re.S).group(1)
 css='\n'.join(x.rstrip()for x in css.splitlines())+'\n.files{display:grid;gap:.65rem}.scroll{overflow-x:auto}td,th{vertical-align:top}.proof-links{display:flex;gap:1rem;flex-wrap:wrap;margin:1rem 0}.math[data-display=false]{font-size:1em}'
 b='<p class="eyebrow">Anticoherent subspaces · Quantum error correction and joint compression</p>'+link('../t2-nesting-upper/index.html','← Earlier nesting upper-bound analysis')
 b+='<h1>Testing stronger second-order upper bounds</h1>'
 b+=p('<strong>No stronger unrestricted upper bound is proved here.</strong> The investigation tests intrinsic quantum-code enumerators and arbitrary matrix-valued compression witnesses. It also gives a rigorous conditional refinement when every polynomial in the subspace has a common factor. The catalogue intervals remain unchanged.')
 b+='<div class="result"><strong>Current general upper bound</strong>'+math(r'K(N,2)\le\min\!\left\{\left\lfloor\frac{N-2}{4}\right\rfloor,\left\lfloor\frac{N+5}{5}\right\rfloor\right\},\qquad N\ge8.')+'</div>'
 b+='<nav class="jump-links" aria-label="On this page"><a href="#qec">Quantum-code relaxation</a><a href="#matrix">Matrix-valued witnesses</a><a href="#factor">Common polynomial factors</a><a href="#files">Proofs and verification</a></nav>'
 b+='<section id="qec"><h2>Intrinsic quantum-code enumerators</h2>'
 b+=p('For a rank-'+math('k',False)+' projector '+math('P',False)+' in spin '+math('j=N/2',False)+', second-order anticoherence means')
 b+=math(r'PJ_aP=0,\qquad PJ_aJ_bP=c\delta_{ab}P,\qquad c=\frac{j(j+1)}3.')
 b+=p('These are Knill–Laflamme error-correction identities with a specified isotropic coefficient matrix. The intrinsic MacWilliams transform for permutation-invariant codes provides necessary linear constraints on the squared spin-multipole coefficients. The transform is expressed in Racah polynomials; see '+link('https://arxiv.org/abs/2605.15372','Teixeira’s 2026 MacWilliams-transform preprint')+' and '+link('https://arxiv.org/html/2605.22439v1#S9','Kubischta–Teixeira, 2026, Section IX')+'.')
 b+=p('Using Hilbert–Schmidt orthonormal spin tensors '+math('T_{\ell m}',False)+', set')
 b+=math(r'q_\ell=\frac1k\sum_{m=-\ell}^{\ell}|\operatorname{Tr}(PT_{\ell m})|^2,\qquad q_\ell\ge0,\quad\sum_{\ell=0}^Nq_\ell=1,\quad q_0=\frac{k}{N+1}.')
 b+=p('The exact rational transform used in the calculation is')
 b+=math(r'M_{ba}=\frac{2b+1}{N+1}\,{}_4F_3\!\left(\begin{matrix}-b,b+1,-a,a+1\\1,-N,N+2\end{matrix};1\right).')
 b+=p('The basic constraints include')
 b+=math(r'q_1=q_2=0,\qquad(Mq)_1=(Mq)_2=0,\qquad Mq\ge0,\qquad q_\ell\le k(Mq)_\ell.')
 b+=p('The tested relaxation also includes time-reversal shadow positivity and capacity constraints. It retains additional transition and shadow inequalities between the mutually orthogonal code, spin-error, and complementary projectors')
 b+=math(r'P,\qquad E=\frac1c\sum_aJ_aPJ_a,\qquad R=I-P-E,\qquad\operatorname{rank}(E)=3k,\quad\operatorname{rank}(R)=N+1-4k.')
 b+=p('The full list and its normalization are given in the '+link('intrinsic-lp-proof.tex','intrinsic-enumerator proof note')+'. Every transform coefficient is computed rationally. The independent verifier reconstructs the transform from factorial six-j sums instead of importing the construction’s hypergeometric recurrence.')
 b+='<h3>Three exact feasible relaxation points</h3><div class="scroll"><table><thead><tr><th>Qubits</th><th>Tested dimension</th><th>Exact conclusion</th><th>Certificate</th></tr></thead><tbody>'
 for r in audit['records']:
  b+=f'<tr><td>{r["n"]}</td><td>{r["k"]}</td><td>All tested necessary inequalities are feasible</td><td>'+link(r['pointFile'],'Rational point')+'</td></tr>'
 b+='</tbody></table></div>'
 b+='<div class="scope-note">'+p('<strong>This is not a physical existence certificate.</strong> The rational points certify that this particular relaxation cannot exclude the tested dimensions. They do not construct a projector or an AC subspace, and they do not prove that these dimensions are attainable. In particular, neither the lower bounds at 35 and 42 qubits nor any catalogue upper bound changes.')+'</div>'
 b+='<div class="proof-links">'+link('verify_lp_feasible.py','Independent exact verifier')+link('lp-feasibility-independent-audit.json','Verification report')+link('intrinsic-lp-scan.json','Exploratory numerical scan')+'</div>'
 b+=p('The wider floating-point scan is exploratory. Its feasible outputs are not promoted to exact certificates beyond the three independently audited points above, and a finite scan gives no all-N feasibility theorem.')+'</section>'
 b+='<section id="matrix"><h2>Arbitrary matrix-valued joint-compression witnesses</h2>'
 b+=p('Let '+math('A_1,\ldots,A_8',False)+' span the traceless spin tensors of ranks one and two. For any ancillary dimension '+math('a',False)+' and arbitrary Hermitian ancillary matrices '+math('H_\mu',False)+', define')
 b+=math(r'W=\sum_{\mu=1}^8A_\mu\otimes H_\mu.')
 b+=p('Every AC2 code satisfies '+math('(P\otimes I_a)W(P\otimes I_a)=0',False)+'. Its dimension is therefore bounded by the inertia of '+math('W',False)+':')
 b+=math(r'\boxed{ak\le n_0(W)+\min\{n_+(W),n_-(W)\}.}')
 b+=p('This family permits general ancillary matrix coefficients, beyond the earlier rotationally invariant scalar quadratics. Alternating-projection searches at 35 qubits attempted to exclude dimension eight:')
 b+='<div class="scroll"><table><thead><tr><th>Ancillary dimension</th><th>Target positive index</th><th>Starts × iterations</th><th>Best boundary eigenvalue</th></tr></thead><tbody>'
 for r in search['experiments']:
  b+=f'<tr><td>{r["auxiliaryDimension"]}</td><td>{r["targetPositiveIndex"]}</td><td>{r["trialCount"]} × {r["stepsPerTrial"]}</td><td>{r["bestBoundaryEigenvalue"]:.7f}</td></tr>'
 b+='</tbody></table></div>'
 b+=p('No strict witness was found. The boundary eigenvalues remained positive. A limit with zero eigenvalues would not suffice without accounting for the kernel term in the displayed bound. These nonconvex search failures do not exclude the existence of better witnesses.')
 b+='<div class="proof-links">'+link('matrix-inertia-framework.tex','Universal witness argument')+link('search_matrix_inertia.py','Search implementation')+link('matrix-inertia-search-report.json','Finite search report')+'</div></section>'
 b+='<section id="factor"><h2>A conditional refinement from common polynomial factors</h2>'
 b+=p('Identify a normalized Dicke coefficient vector with the homogeneous binary form')
 b+=math(r'f(x,y)=\sum_{r=0}^Nc_r\sqrt{\binom Nr}\,x^{N-r}y^r.')
 b+=p('For any second-order AC subspace '+math('\mathcal C',False)+' and '+math('N\ge4',False)+', the highest-spin projection gives an injective multiplication map')
 b+=math(r'\mathcal C\otimes\operatorname{Sym}^4\mathbb C^2\longrightarrow\operatorname{Sym}^{N+4}\mathbb C^2,\qquad f\otimes g\longmapsto fg.')
 b+=p('Consequently, if all forms in '+math('\mathcal C',False)+' share a polynomial factor of degree '+math('g',False)+', the image lies in a space of dimension '+math('N+5-g',False)+'. This proves')
 b+=math(r'\boxed{5k\le N+5-g\qquad\text{when the subspace has a common factor of degree }g.}')
 b+=p('Here the common factor is shared by every vector in the subspace, not merely a repeated Majorana point of one chosen representative. In general '+math('g=0',False)+' is possible, so this is a conditional strengthening and does not improve the unrestricted bound.')
 b+=p('For odd '+math('N',False)+', a subspace contained entirely in one Dicke-parity sector has a common linear factor. Hence '+math('5k\le N+4',False)+'. At '+math('N=35',False)+', this excludes dimension eight within either such sector. It does not cover a code using both parity sectors or an arbitrary subspace.')
 b+='<details><summary>Why multiplication by degree-four forms is injective</summary>'
 b+=p('With an ancillary spin two, the established witness')
 b+=math(r'W=X^2+(j+3)X-2j(j+1)I,\qquad X=\mathbf J\cdot\mathbf S')
 b+=p('has zero compression on '+math('\mathcal C\otimes V_2',False)+' and is strictly negative on every total-spin sector below '+math('j+2',False)+'. Thus no nonzero vector of that product space can lie in the kernel of the highest-spin projection. The highest-spin map is the Cartan product, which is binary-form multiplication up to nonzero normalization. The existing '+link('../t2-spin2-inertia/proof.tex','spin-two spectrum proof')+' supplies the all-N signs.')
 b+=p('Equivalently, no nontrivial relation '+math('\sum_af_ag_a=0',False)+' can have homogeneous coefficient forms of common degree at most four. No corresponding degree-five injectivity result is proved; it would imply the stronger bound '+math('6k\le N+6',False)+'.')+'</details>'
 b+='<div class="proof-links">'+link('hankel-syzygy-consequence.tex','Common-factor theorem and proof')+link('../t2-spin2-inertia/verify.py','Existing exact spectral verifier')+'</div></section>'
 b+='<section id="files"><h2>Proofs, scope and reproducibility</h2><div class="files">'
 for path,label in [('intrinsic-lp-proof.tex','Intrinsic quantum-code relaxation and normalization'),('verify_lp_feasible.py','Independent rational verification of the three feasible points'),('certify_lp_feasible.py','Rational feasible-point reconstruction'),('intrinsic_lp.py','Intrinsic-enumerator scan implementation'),('matrix-inertia-framework.tex','General matrix-valued inertia theorem'),('hankel-syzygy-consequence.tex','Conditional common-factor refinement'),('https://arxiv.org/abs/2605.15372','Teixeira: intrinsic MacWilliams transform, 2026'),('https://arxiv.org/abs/2605.22439','Kubischta–Teixeira: permutation-invariant code LP framework, 2026')]:b+=link(path,label)
 b+='</div>'+p('The exact results here are a limitation of the specified linear relaxation at three candidate parameters, the general inertia witness framework, and the conditional common-factor theorem. They do not prove optimality of the current unrestricted bound.')+'</section><footer>'+link('../../index.html#t=2','Return to the AC-subspace catalogue')+'<span>No catalogue dimensions or interval endpoints changed.</span></footer>'
 head='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>QEC and joint-compression upper bounds</title><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script defer src="../../vendor/katex/katex.min.js"></script><style>'+css+'</style><script>document.addEventListener("DOMContentLoaded",()=>document.querySelectorAll("[data-tex]").forEach(el=>{if(window.katex)el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==="true",throwOnError:false,strict:"ignore",trust:false});}));</script></head><body><main>'
 (HERE/'index.html').write_text(head+b+'</main></body></html>\n')
if __name__=='__main__':main()
