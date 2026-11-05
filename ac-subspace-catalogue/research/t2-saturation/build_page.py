#!/usr/bin/env python3
"""Offline research note on restricted cyclic saturation of AC2 bounds."""
from pathlib import Path
import html,json,re
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def tex(s):return '<div class="math" data-tex="'+html.escape(s,quote=True)+'" data-display="true"></div>'
def p(s):return '<p>'+s+'</p>'
def main():
 cert=json.loads((HERE/'construction-cyclic-exclusion.json').read_text())
 assert not cert['general_subspace_exclusion']
 css=re.search(r'<style>(.*?)</style>',(ROOT/'research/t2-spin2-inertia/index.html').read_text(),re.S).group(1)
 css='\n'.join(line.rstrip() for line in css.split('\n'))
 css+=' .math{overflow-x:auto;max-width:100%;padding:.5rem 0}nav{display:flex;gap:1rem;flex-wrap:wrap;margin:1.5rem 0}.scope{border-left:4px solid #8b6339;padding:.25rem 1.1rem;margin:1.5rem 0;background:#f8f3e7}.scroll-table{overflow-x:auto}table{width:100%;border-collapse:collapse}td,th{text-align:left;border-bottom:1px solid #deddd3;padding:.6rem}details{margin:1.4rem 0}summary{cursor:pointer;font-weight:600}.downloads{display:grid;grid-template-columns:1fr 1fr;gap:.6rem}@media(max-width:600px){.downloads{grid-template-columns:1fr}}'
 body='<p class="eyebrow">Anticoherent subspaces · An exclusion for cyclic-invariant constructions</p><a href="../t2-general/index.html">← General construction prescriptions</a><h1>Can the cyclic prescription attain the bound indefinitely?</h1>'
 body+=p('A subspace invariant under a rotation through 2π/k cannot attain dimension k = floor(N/5) + 1 from N = 420 onward. An exact sector argument and projector inequality establish this for any choice of rotation axis.')
 body+='<div class="scope">'+p('<strong>Scope of the result:</strong> this excludes every Cₖ-invariant subspace at the stated target dimension, about any axis. Arbitrary subspaces without that symmetry remain unexcluded. No infinite family of general subspaces attaining the bound has been established here.')+'</div>'
 body+='<nav aria-label="On this page"><a href="#theorem">Exact theorem</a><a href="#sectors">Why every residue has rank one</a><a href="#projector">Projector inequality</a><a href="#window">Seven-residue window</a><a href="#contradiction">Moment contradiction</a><a href="#subsequence">Sharper subsequence</a><a href="#searches">Numerical exploration</a><a href="#files">Proof files</a></nav>'
 body+='<section id="theorem"><h2>The symmetry hypothesis is essential</h2>'+p('Let')+tex(r'N\ge420,\qquad k=\left\lfloor\frac N5\right\rfloor+1.')
 body+=p('There is no k-dimensional second-order AC subspace invariant under a rotation through the angle below, about any axis:')+tex(r'\phi=\frac{2\pi}{k}.')
 body+=p('Rotate that axis to z. Invariance decomposes the code into Dicke residue eigenspaces modulo k. The sector lemma below forces exactly one code dimension in every residue, so the code must admit a basis of the form')+tex(r'|u_{a+1}\rangle=\sum_{\substack{0\le r\le N\\r\equiv a\pmod k}}c_r|D_N^{(r)}\rangle,\qquad a=0,\ldots,k-1.')
 body+=p('Coefficients may be complex, may vanish, and may vary arbitrarily between neighboring residues. No smoothness or phase-continuity assumption is made. The rank-one residue conclusion is proved, rather than assumed.')
 body+='<div class="result"><strong>Cyclic-invariant exclusion</strong>'+tex(r'\boxed{\text{No }C_k\text{-invariant AC2 subspace has dimension }k=\lfloor N/5\rfloor+1\text{ for }N\ge420.}')+'</div>'
 body+=p('The exact five-dimensional example at N = 24 remains valid and optimal. The theorem prevents an infinite sequence of such optimal examples with Cₖ symmetry at the target dimension. <a href="../t2-general/n24-k5/index.html">Read the N = 24 construction ↗</a>')+'</section>'
 body+='<section id="sectors"><h2>Why every residue sector must have rank one</h2>'+p('Write the target relation and isotropic variance as')+tex(r'N=5k-s,\qquad1\le s\le5,\qquad c=\frac{N(N+2)}{12}.')
 body+=p('Each residue contains at most five equally spaced spin coordinates. Pad a shorter list if necessary and write them as')+tex(r'd-2k,\quad d-k,\quad d,\quad d+k,\quad d+2k.')
 body+=p('Suppose the AC subspace inside one sector had dimension at least two. There is then a nonzero vector whose top coefficient vanishes, and another whose bottom coefficient vanishes. Both have mean zero and variance c. A mean-zero distribution supported in an interval has variance at most minus the product of its endpoints, so these two vectors give')+tex(r'\begin{aligned}c&\le-(d-2k)(d+k)=2k^2+kd-d^2,\\c&\le-(d-k)(d+2k)=2k^2-kd-d^2.\end{aligned}')
 body+=p('Consequently any sector of dimension at least two would require')+tex(r'c\le2k^2-|kd|-d^2\le2k^2.')
 body+=p('But in the theorem’s range k is at least 85, and')+tex(r'12(c-2k^2)=k^2+10(1-s)k+s(s-2)\ge k^2-40k+15>0.')
 body+=p('Thus each of the k residue sectors has dimension at most one. Since the full code has dimension k, every residue must contribute exactly one vector. This reduces all Cₖ-invariant target-dimensional codes to the support structure used below.')+'</section>'
 body+='<section id="projector"><h2>Four orthogonal error images</h2>'+p('Write')+tex(r'j=N/2,\qquad c=\frac{j(j+1)}3,\qquad P=\text{projector onto the AC subspace}.')
 body+=p('Second-order anticoherence gives')+tex(r'PJ_\alpha P=0,\qquad PJ_\alpha J_\beta P=c\delta_{\alpha\beta}P.')
 body+=p('The identity and the three spin operators, normalized on the code, therefore have mutually orthogonal images. Their combined range projector satisfies')+tex(r'E=P+\frac1c\sum_{\alpha=x,y,z}J_\alpha PJ_\alpha\le I.')
 body+=p('Set the Dicke-basis diagonal entries and spin coordinates to')+tex(r'p_r=\langle D_N^{(r)}|P|D_N^{(r)}\rangle,\qquad m_r=j-r.')
 body+=p('Taking a diagonal entry of the projector inequality yields')+tex(r'\left(1+\frac{m_r^2}{c}\right)p_r+\frac{r(N-r+1)p_{r-1}+(r+1)(N-r)p_{r+1}}{2c}\le1.\tag{1}')
 body+=p('Terms outside the physical index range are zero. This inequality holds for every AC2 subspace. The rank-one sector conclusion is used in the next step: every residue separately has unit weight and the fixed AC moments.')+tex(r'\sum_{r\equiv a\ (k)}p_r=1,\qquad\sum_{r\equiv a\ (k)}m_rp_r=0,\qquad\sum_{r\equiv a\ (k)}m_r^2p_r=c.\tag{2}')+'</section>'
 body+='<section id="window"><h2>A seven-residue window bounds each band weight</h2>'+p('Average the first seven residues with positive sine weights:')+tex(r'\theta=\frac\pi8,\qquad f_a=\sin((a+1)\theta),\quad0\le a\le6,\qquad F=\sum_{a=0}^6f_a,\quad f_{-1}=f_7=0.')
 body+=p('Write')+tex(r'N=5k-s,\qquad1\le s\le5.')
 body+=p('All five bands contain these seven indices because the stated range guarantees k at least s + 6. In band ℓ, multiply inequality (1) at r = ℓk + a by fₐ and sum over the window. The coefficient of an internal diagonal entry is')
 body+=tex(r'\begin{aligned}&\left(1+\frac{m_r^2}{c}\right)f_a+\frac{(j(j+1)-m_r^2+m_r)f_{a+1}+(j(j+1)-m_r^2-m_r)f_{a-1}}{2c}\\&=\left[1+3\cos\theta+\frac{m_r^2}{c}(1-\cos\theta)\right]f_a+\frac{m_r}{c}\sin\theta\cos((a+1)\theta)\\&\ge\left[1+3\cos\theta-\frac6{N+2}\right]f_a.\end{aligned}')
 body+=p('The last step follows from')+tex(r'\sin\theta\,|\cot((a+1)\theta)|\le1,\qquad\frac{|m_r|}{c}\le\frac6{N+2}.')
 body+=p('All terms from outside the window are nonnegative. Use the exact lower bound')+tex(r'\cos(\pi/8)>\frac{11}{12}.')
 body+=p('It follows that the five averaged band weights satisfy')+tex(r'\begin{gathered}B(N)=\frac{15}{4}-\frac6{N+2},\qquad W_\ell=\frac1F\sum_{a=0}^6f_a p_{\ell k+a},\\0\le W_\ell\le\frac1{B(N)},\qquad\sum_{\ell=0}^4W_\ell=1.\tag{3}\end{gathered}')
 body+=p('This step controls even rapidly oscillating coefficients. It does not replace them by a smooth limiting profile.')+'</section>'
 body+='<section id="contradiction"><h2>A quadratic moment cannot fit those band caps</h2>'+p('Use the five ideal normalized latitudes and the quadratic polynomial')+tex(r'z_\ell=1-\frac{2\ell}{5},\qquad Q(z)=z^2-\frac85z.')
 body+=p('The actual normalized spin coordinate differs by a controlled nonnegative amount:')+tex(r'\frac{m_{\ell k+a}}j=z_\ell-e_{\ell a},\qquad e_{\ell a}=\frac{2(a+s\ell/5)}N\le\frac{12+2\ell}{N}=:\epsilon_\ell.')
 body+=p('At the ideal latitudes, the excess above minus three fifths is')+tex(r'\left(Q(z_\ell)+\frac35\right)_{\ell=0}^4=\left(0,0,\frac8{25},\frac{24}{25},\frac{48}{25}\right).')
 body+=p('In the first band the correction does not raise the polynomial above minus three fifths. In each remaining band,')+tex(r'Q(z_\ell-e)\le Q(z_\ell)+\left(\frac85-2z_\ell\right)\epsilon_\ell+\epsilon_\ell^2,\qquad1\le\ell\le4.')
 body+=p('The four excesses are nonnegative, so the band cap (3) gives the following upper bound for the averaged polynomial moment:')+tex(r'\begin{gathered}\langle Q\rangle=\frac1F\sum_{a=0}^6f_a\sum_{\ell=0}^4p_{\ell k+a}\,Q\!\left(\frac{m_{\ell k+a}}j\right),\\\langle Q\rangle\le U(N):=-\frac35+\frac{\displaystyle\frac{16}{5}+\frac{584}{5N}+\frac{1176}{N^2}}{\displaystyle\frac{15}{4}-\frac6{N+2}}.\tag{4}\end{gathered}')
 body+=p('But the exact zero first moment and isotropic second moment in every residue, equation (2), give')+tex(r'\langle Q\rangle=\frac{c}{j^2}=\frac{N+2}{3N}>\frac13.\tag{5}')
 body+=p('The positive numerator in (4) decreases with N, while its positive denominator increases. Consequently U decreases throughout the stated range. At its endpoint, exact rational arithmetic gives')+tex(r'\frac13-U(420)=\frac{872}{1655325}>0.')
 body+=p('Thus (4) requires the moment to be below one third, while (5) requires it to be above one third. This proves the cyclic-invariant exclusion for every N at least 420.')+'</section>'
 body+='<section id="subsequence"><h2>A sharper threshold when N = 5k − 1</h2>'+p('Here s = 1, so the latitude correction improves to')+tex(r'\epsilon_\ell=\frac{12+2\ell/5}{N}.')
 body+=p('The upper moment bound becomes')+tex(r'U_1(N)=-\frac35+\frac{\displaystyle\frac{16}{5}+\frac{424}{5N}+\frac{3384}{5N^2}}{\displaystyle\frac{15}{4}-\frac6{N+2}}.')
 body+=p('At the first threshold in this subsequence,')+tex(r'N=309=5\cdot62-1,\qquad\frac13-U_1(309)=\frac{93358}{738545535}>0.')
 body+=p('On this subsequence, c − 2k² = (k² − 1)/12 is also positive, so the sector lemma applies. The same monotonicity argument excludes every Cₖ-invariant dimension-k AC2 subspace for k at least 62, about any axis. Arbitrary subspaces without this symmetry remain unexcluded.')+'</section>'
 body+='<section id="searches"><h2>Numerical exploration is separate from the proof</h2>'+p('Before deriving the obstruction, searches tested cyclic and reflection-symmetric extensions of the exact N = 24 example. Random starts and interpolation of that basis did not produce an exact candidate at the following larger target dimensions.')
 body+='<details><summary>Search residuals and their limitations</summary><div class="scroll-table"><table><thead><tr><th>N</th><th>Target k</th><th>Best normalized residual</th></tr></thead><tbody>'
 for n,k,res in [(29,6,'3.69 × 10⁻⁴'),(34,7,'2.69 × 10⁻²'),(39,8,'3.01 × 10⁻²'),(44,9,'4.13 × 10⁻²'),(49,10,'5.46 × 10⁻²')]:body+=f'<tr><td>{n}</td><td>{k}</td><td>{res}</td></tr>'
 body+='</tbody></table></div>'+p('These are local optimization outcomes, not lower bounds on the attainable residual. They prove no exclusion, even within the ansatz. In particular, the theorem above does not settle these smaller qubit numbers.')+p('<a href="cyclic-extension-search-report.json">Search methods, counts, and limitations ↓</a>')+'</details>'
 body+=p('The exact exclusion uses only the projector inequality, residue moments, the finite sine window, and rational endpoint gaps. It does not depend on the search residuals.')+'</section>'
 body+='<section id="files"><h2>Proof and verification</h2><div class="downloads">'
 files=[('construction-cyclic-exclusion.tex','Complete LaTeX proof'),('construction_cyclic_exclusion.py','Exact arithmetic verifier'),('construction-cyclic-exclusion.json','Verified coefficients and thresholds'),('cyclic-exclusion-independent-audit.py','Independent exact verifier'),('cyclic-exclusion-independent-audit.json','Independent audit report'),('../t2-general/index.html','General analytical lower-bound constructions')]
 for file,label in files:body+='<a href="'+file+'">'+label+' ↗</a>'
 body+='</div>'+p('The unrestricted problem remains open in this work: a subspace without the specified Cₖ symmetry could evade the obstruction. The established general construction families continue to provide rigorous lower bounds.')+'</section><footer><p>Anticoherent subspaces · Research note</p><a href="../../index.html">Return to the catalogue ↗</a></footer>'
 head='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>An exclusion for cyclic-invariant saturation of AC-subspace bounds</title><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script defer src="../../vendor/katex/katex.min.js"></script><style>'+css+'</style><script>document.addEventListener("DOMContentLoaded",()=>document.querySelectorAll("[data-tex]").forEach(el=>{if(window.katex)el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:true,throwOnError:false,strict:"ignore",trust:false});}));</script></head><body><main>'
 (HERE/'index.html').write_text(head+body+'</main></body></html>\n')
if __name__=='__main__':main()
