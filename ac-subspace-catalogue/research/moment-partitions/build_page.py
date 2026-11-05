#!/usr/bin/env python3
"""Build the prescribed-moment construction research note from exact records."""
from pathlib import Path
import html,json,re
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def math(x,display=True):
    tag='div'if display else'span'
    return f'<{tag} class="math" data-display="{str(display).lower()}" data-tex="{html.escape(x,quote=True)}"></{tag}>'
def p(x):return '<p>'+x+'</p>'
def link(path,label):return f'<a href="{path}">{label} ↗</a>'
def links(items):return '<div class="proof-links">'+''.join(link(a,b)for a,b in items)+'</div>'
def main():
    family=json.loads((HERE/'families-fourband-certificate.json').read_text());assert family['passed']
    independent=json.loads((HERE/'search-independent-family-audit.json').read_text());assert independent['passed']
    examples=json.loads((HERE/'search-exact-examples.json').read_text())['subspaces']
    assert all(x['verification']['allCrossReductionsVerified']for x in examples)
    geom=json.loads((HERE/'geometry_verification.json').read_text());assert all(x['passed']for x in geom)
    css=re.search(r'<style>(.*?)</style>',(ROOT/'research/t2-qec-upper/index.html').read_text(),re.S).group(1)
    css+='\nfigure{margin:2rem 0}figure img{width:100%;height:auto;border:1px solid var(--line,#d7d9d3);border-radius:4px}figcaption{font-size:.92rem;color:var(--muted,#5c6664);line-height:1.6;margin-top:.7rem}.math{max-width:100%;overflow-x:auto}.formula-set{padding:.8rem 0}.keyline{font-weight:600}.steps{line-height:1.8}.steps li{margin:.7rem 0}.note{font-size:.92rem;color:var(--muted,#5c6664)}'
    b='<p class="eyebrow">Anticoherent subspaces · Prescribed moments and analytical constructions</p>'
    b+=link('../../index.html#t=2','← AC-subspace catalogue')+'<h1>Constructing AC2 subspaces from prescribed moments</h1>'
    b+=p('Convex combinations of lifted Dicke levels supply exact first and second moments. Adding controlled phase cancellation gives a stronger general construction, with every coefficient specified by rational numbers and finite radical recurrences.')
    b+='<div class="result"><strong>Improved constructive lower bound</strong>'+math(r'\boxed{K(N,2)\ge\left\lfloor\frac{N-1}{8}\right\rfloor,\qquad N\ge17.}')+p('The proof covers every integer in this range. It increases the earlier four-band construction by one dimension whenever '+math(r'N\not\equiv0\pmod8',False)+'. It does not establish that these subspaces have maximum possible dimension.')+'</div>'
    b+='<nav class="jump-links" aria-label="On this page"><a href="#moments">Moment geometry</a><a href="#family">General radical family</a><a href="#examples">Exact examples</a><a href="#overlap">Overlapping supports</a><a href="#limits">Sparse-support limit</a><a href="#files">Proofs and verification</a></nav>'
    b+='<section id="moments"><h2>The target is a prescribed point on the moment plane</h2>'
    b+=p('Put '+math(r'j=N/2,\ m=j-r,\ c=j(j+1)/3',False)+'. A probability distribution on Dicke indices has the required diagonal AC2 moments when')
    b+=math(r'\sum_r w_r=1,\qquad \sum_r w_r m_r=0,\qquad \sum_r w_r m_r^2=c.')
    b+=p('Equivalently, the lifted points '+math(r'(m_r,m_r^2)',False)+' have a convex combination equal to '+math(r'(0,c)',False)+'. If the used Dicke indices are globally separated by at least three, all spin operators shifting excitation by one or two vanish on their coordinate space. Disjoint convex blocks then give an orthonormal AC2 basis:')
    b+=math(r'|u_a\rangle=\sum_{r\in B_a}\sqrt{w_r}\,|D_N^{(r)}\rangle.')
    b+=p('For three distinct magnetic levels '+math(r'a,b,d',False)+', the probabilities are explicit:')
    b+=math(r'w_a=\frac{c+bd}{(a-b)(a-d)},\quad w_b=\frac{c+ad}{(b-a)(b-d)},\quad w_d=\frac{c+ab}{(d-a)(d-b)}.')
    b+=p('A two-point block exists exactly when '+math(r'ab=-c',False)+'. The weights are '+math(r'-b/(a-b)',False)+' and '+math(r'a/(a-b)',False)+'. Positivity is the geometric condition that the target lies in the block’s segment or triangle.')
    b+='<figure><img src="moment-partition-n22.svg" alt="Two segments and one triangle of lifted magnetic levels containing the common target (0,44) for a 22-qubit AC2 basis."><figcaption>At 22 qubits, two pairs and one triangle all contain the same prescribed moment point. Their disjoint Dicke supports yield the real radical basis below.</figcaption></figure>'
    b+='<details><summary>The exact number of disjoint blocks on a fixed sparse support</summary>'
    b+=p('Let '+math(r'M',False)+' be the number of available nodes, '+math(r'd',False)+' the Tukey depth of the target, and '+math(r'p',False)+' the number of pairs satisfying '+math(r'ab=-c',False)+'. Those pairs are automatically disjoint. The largest number of disjoint moment blocks is')
    b+=math(r'\min\!\left\{d,\left\lfloor\frac{M+p}{3}\right\rfloor\right\}.')
    b+=p('The depth counts the fewest nodes in a generic halfplane through the target. Carathéodory’s theorem reduces blocks to pairs or triangles; the prescribed planar Tverberg theorem supplies a matching construction after the pairs are removed. This statement concerns a fixed support set and disjoint blocks, not arbitrary AC2 subspaces. See the '+link('geometry.tex','geometry proof')+' and '+link('https://doi.org/10.4230/LIPIcs.ESA.2021.51','Har-Peled–Zhou, Theorem 7')+'.')+'</details></section>'
    b+='<section id="family"><h2>Four Dicke bands with phases: a stronger general family</h2>'
    b+=p('A globally three-separated support leaves too few available Dicke levels for a faster construction. The new family permits levels two apart inside each of four bands. Convex moment weights solve the diagonal conditions, and four-sided polygons provide the phases that cancel every remaining cross term.')
    b+=math(r'q=\left\lfloor\frac{N-1}{8}\right\rfloor,\qquad N=8q+\rho,\quad1\le\rho\le8,\qquad k=q.')
    b+=p('For '+math(r'1\le\rho\le7',False)+', take '+math(r'h=2q+1,\ s=\lfloor(\rho-1)/2\rfloor',False)+'. For '+math(r'\rho=8',False)+', take '+math(r'h=2q+2,\ s=2',False)+'. The basis has the explicit form')
    b+=math(r'|u_i\rangle=\sum_{\ell=0}^{3}\sqrt{w_{i\ell}}\,\eta_{i\ell}|D_N^{(s+2i+\ell h)}\rangle,\qquad 0\le i<q.')
    b+=p('The weights are positive rational numbers. The unit phases '+math(r'\eta_{i\ell}',False)+' are obtained recursively by drawing a quadrilateral of prescribed side lengths; its diagonal and angles require only field operations and square roots. No optimization or numerical root fitting is part of this prescription.')
    b+='<details><summary>Rational weights and exact radical phase recurrence</summary><div class="formula-set">'
    b+=math(r'x_i=\frac{j-s-2i}{h},\qquad v=\frac{j(j+1)}{3h^2},\qquad \epsilon=\begin{cases}1/h,&\rho<8,\\0,&\rho=8.\end{cases}')
    b+=math(r'\begin{aligned}b_0(x)&=\frac{v+x^2}{4}-\frac{35x}{36}+\frac56,&b_1(x)&=-\frac{v+x^2}{4}+\frac{5x}{12}+\frac12,\\b_2(x)&=-\frac{v+x^2}{4}+\frac{13x}{12}-\frac12,&b_3(x)&=\frac{v+x^2}{4}-\frac{19x}{36}+\frac16.\end{aligned}')
    b+=math(r'(w_{i0},w_{i1},w_{i2},w_{i3})=(b_0(x_i),b_1(x_i),b_2(x_i),b_3(x_i))+\epsilon(x_i-3/2)(-1/3,1,-1,1/3).')
    b+=p('For adjacent rows define '+math(r'r_{i\ell}=s+2i+\ell h',False)+' and')
    b+=math(r'L_{i\ell}=\sqrt{w_{i\ell}w_{i+1,\ell}(r_{i\ell}+1)(r_{i\ell}+2)(N-r_{i\ell})(N-r_{i\ell}-1)}.')
    b+=p('These four positive lengths satisfy strict quadrilateral inequalities for every allowed '+math(r'N,i',False)+'. Write them as '+math(r'(a,b,c,d)',False)+' and choose')
    b+=math(r'D=\frac{\max\{|a-b|,|c-d|\}+\min\{a+b,c+d\}}2.')
    b+=math(r'Z(u,v)=\left(C_u+\mathrm i\sqrt{1-C_u^2},C_v-\mathrm i\sqrt{1-C_v^2}\right),\quad C_u=\frac{D^2+u^2-v^2}{2Du},\quad C_v=\frac{D^2+v^2-u^2}{2Dv}.')
    b+=math(r'(z_0,z_1)=Z(a,b),\qquad(z_2,z_3)=-Z(c,d),\qquad \eta_{0\ell}=1,\quad \eta_{i+1,\ell}=\eta_{i\ell}z_\ell.')
    b+=p('Then '+math(r'\sum_\ell L_{i\ell}z_\ell=0',False)+'. The bands have gaps three or four; the only distance-two cross entries connect adjacent rows within a band. The displayed recurrence cancels all of them.')+'</div></details>'
    b+='<h3>Improvement over the preceding general formulas</h3>'
    b+=p('The table compares the strongest preceding general construction '+math(r'\max\{\lfloor(N+3)/9\rfloor,\lfloor N/8\rfloor-1\}',False)+' with the new family. Special constructions can still give larger dimensions at particular sizes.')
    b+='<div class="scroll"><table><thead><tr><th>Qubits '+math(r'N',False)+'</th><th>Previous general lower bound</th><th>New radical family</th></tr></thead><tbody>'
    for n in [41,49,57,65]:b+=f'<tr><td>{n}</td><td>{max((n+3)//9,n//8-1)}</td><td>{(n-1)//8}</td></tr>'
    b+='</tbody></table></div>'
    b+=p('Combining the earlier separated triples with the new theorem gives, for every '+math(r'N\ge17',False)+',')
    b+=math(r'K(N,2)\ge\max\!\left\{\left\lfloor\frac{N+3}{9}\right\rfloor,\left\lfloor\frac{N-1}{8}\right\rfloor\right\}.')
    b+=p('Exact polynomial coefficient certificates establish positivity and phase closure for all eight remainders. An independent implementation reconstructs the weights from the physical moment equations and verifies the same all-integer domain. Direct spin-matrix checks for 144 sizes '+math(r'17\le N\le160',False)+' separately validate the numerical rendering.')
    b+=links([('families-fourband-proof.tex','All-N theorem and proof'),('families-fourband-certificate.json','Exact coefficient certificate'),('search-independent-family-audit.json','Independent all-N audit'),('families-example-weights.json','Exact example data, including N=41, k=5')])+'</section>'
    b+='<section id="examples"><h2>Two particularly simple analytical bases</h2>'
    b+=p('These bases simplify or illustrate already attained dimensions. They do not change the catalogue’s unrestricted lower or upper endpoints.')
    for record in examples:
        n=record['n'];k=record['dimension'];ident=record['id']
        b+=f'<h3>{n} qubits · dimension {k}</h3>'
        if n==12:
            b+=p('Five globally separated levels support a two-dimensional real radical compression. The two vectors overlap on the same support. No disjoint partition into two moment blocks exists on these five nodes, so this example shows what overlap adds.')
            formulas=[r'|u\rangle=\sqrt{\frac{77}{486}}\bigl(|D_{12}^{(0)}\rangle+|D_{12}^{(12)}\rangle\bigr)+\sqrt{\frac{35}{243}}\bigl(|D_{12}^{(3)}\rangle+|D_{12}^{(9)}\rangle\bigr)+\frac{4\sqrt2}{9}|D_{12}^{(6)}\rangle.',r'|w\rangle=\sqrt{\frac5{54}}\bigl(|D_{12}^{(0)}\rangle-|D_{12}^{(12)}\rangle\bigr)-\sqrt{\frac{11}{27}}\bigl(|D_{12}^{(3)}\rangle-|D_{12}^{(9)}\rangle\bigr).']
            b+=math(r'S=\{6,3,0,-3,-6\},\qquad c=14,\qquad K(12,2)=2\ \text{(already proved)}.')
        else:
            b+=p('Two pairs and one triangle give a three-dimensional basis with positive real radical coefficients. The seven used Dicke levels have global spacing at least three. This simplifies the earlier complex-phase witness.')
            b+=math(r'B_1=\{0,15\},\qquad B_2=\{7,22\},\qquad B_3=\{3,10,18\},\qquad c=44.')
            formulas=[r'|u_1\rangle=\sqrt{\frac4{15}}|D_{22}^{(0)}\rangle+\sqrt{\frac{11}{15}}|D_{22}^{(15)}\rangle.',r'|u_2\rangle=\sqrt{\frac{11}{15}}|D_{22}^{(7)}\rangle+\sqrt{\frac4{15}}|D_{22}^{(22)}\rangle.',r'|u_3\rangle=\sqrt{\frac{37}{105}}|D_{22}^{(3)}\rangle+\sqrt{\frac3{14}}|D_{22}^{(10)}\rangle+\sqrt{\frac{13}{30}}|D_{22}^{(18)}\rangle.']
        b+='<details><summary>Show every analytical basis coefficient</summary>'+''.join(math(x)for x in formulas)+'</details>'
        b+=p(str(record['verification']['scalarEntriesChecked'])+' ordered two-qubit cross-reduction entries are verified exactly, including every pair of different basis vectors.')
        b+=links([(f'../../index.html#n={n}&t=2&subspace={ident}','Open this representative in the catalogue')])
    b+=links([('search-examples-proof.tex','Both basis derivations'),('search-exact-examples.json','Complete analytical and numerical basis data'),('search-exact-verification.json','Exact cross-reduction verification')])+'</section>'
    b+='<section id="overlap"><h2>An exact construction algorithm allowing overlapping supports</h2>'
    b+=p('On any globally three-separated support, the full AC2 problem reduces to two diagonal Hermitian forms:')
    b+=math(r'Q^*Q=I_k,\qquad Q^*DQ=0,\qquad Q^*(D^2-cI)Q=0.')
    b+=p('Choose a small rational '+math(r'\varepsilon',False)+' so that '+math(r'a_m=m+\varepsilon(m^2-c)',False)+' never vanishes. Sort the ratios '+math(r'(m^2-c)/a_m',False)+' separately on the positive and negative sectors, obtaining '+math(r'\alpha_1<\cdots<\alpha_u',False)+' and '+math(r'\beta_1<\cdots<\beta_v',False)+'. A dimension-'+math(r'k',False)+' compression exists precisely when '+math(r'k\le\min(u,v)',False)+' and')
    b+=math(r'\max(\alpha_i,\beta_i)\le\min(\alpha_{u-k+i},\beta_{v-k+i}),\qquad 1\le i\le k.')
    b+=p('The '+link('https://doi.org/10.4153/CJM-1957-036-1','Fan–Pall interlacing theorem')+' supplies compatible compressed spectra. Successive codimension-one compressions then construct an exact basis using rational arithmetic and square roots. The supplied implementation includes examples at 12 and 21 qubits. This criterion solves the existence question on the chosen support, not the unrestricted AC2 problem.')
    b+=links([('geometry.tex','Moment geometry and graph-compression proof'),('geometry_graph.py','Constructive radical algorithm'),('geometry_examples.json','Algorithmic examples'),('geometry_verification.json','Independent exact cross-reduction checks')])+'</section>'
    b+='<section id="limits"><h2>A rigorous limit of globally sparse supports</h2>'
    b+=p('Even if columns overlap and the available Dicke levels are chosen nonuniformly, global spacing at least three imposes')
    b+=math(r'\boxed{k\le\left\lceil\frac N9\right\rceil\quad\text{for any AC2 code on a globally three-separated support}.}')
    b+=p('A quadratic with zero code compression can be chosen to be positive only on support levels '+math(r'm\ge(j+1)/3',False)+'. Projection into these positive coordinates is injective. The interval has length '+math(r'(N-1)/3',False)+', so at spacing three it contains at most '+math(r'\lceil N/9\rceil',False)+' support points.')
    b+=p('This is a bound on the sparse-support mechanism. The four-band family escapes it by using distance-two levels and canceling the corresponding cross terms with phases. Convex diagonal moments alone would not establish anticoherence in that denser setting.')
    b+=links([('geometry.tex','Sparse-support bound and proof')])+'</section>'
    b+='<section id="files"><h2>Proofs and reproducible data</h2><div class="files">'
    for path,label in [('families-fourband-proof.tex','General radical family: full proof'),('families_fourband.py','All-N coefficient certificate generator'),('search_independent_family_audit.py','Independent all-N physical-moment audit'),('families_render.py','Explicit weights, radical recurrence and numerical rendering'),('families-render-checks.json','Full spin-compression rendering checks'),('search_exact_examples.py','Exact 12- and 22-qubit basis construction and verification'),('geometry.tex','Convex partitions, overlapping compressions and sparse-support limit'),('geometry_graph.py','Fan–Pall construction implementation'),('geometry_verify.py','Independent Dicke-splitting verification')]:b+=link(path,label)
    b+='</div>'+p('The general lower bound is proved for all allowed integers. Exact example certificates check entire subspaces rather than only the individual basis vectors. No unrestricted optimality claim follows from these constructions.')+'</section><footer>'+link('../../index.html#t=2','Return to the AC-subspace catalogue')+'<span>Exact analytical constructions · upper bounds unchanged</span></footer>'
    head='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Prescribed-moment constructions of AC2 subspaces</title><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script defer src="../../vendor/katex/katex.min.js"></script><style>'+css+'</style><script>document.addEventListener("DOMContentLoaded",()=>document.querySelectorAll("[data-tex]").forEach(el=>{if(window.katex)el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==="true",throwOnError:false,strict:"ignore",trust:false});}));</script></head><body><main>'
    (HERE/'index.html').write_text(head+b+'</main></body></html>\n')
if __name__=='__main__':main()
