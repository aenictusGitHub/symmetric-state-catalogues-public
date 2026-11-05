#!/usr/bin/env python3
"""Render the Molien construction notes and changes from the verified exports."""
from pathlib import Path
import html,json,hashlib
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
read=lambda p:json.loads(p.read_text())
data=read(ROOT/'data/catalogue.json')
source=read(HERE/'bases.json')
before={(e['n'],e['t']):e for e in read(HERE/'baseline-entries.json')['entries']}
records={r['id']:r for r in source['subspaces']}
selected=[r for r in data['subspaces'] if r['id'] in records]
changes=[]
for e in data['entries']:
 old=before[e['n'],e['t']]
 assert e['upper']==old['upper'],e['id']
 if e['lower']>old['lower']:
  ids=[r['id'] for r in selected if r['n']==e['n'] and r['t']>=e['t'] and r['dimension']==e['lower']]
  assert ids
  changes.append({'n':e['n'],'t':e['t'],'previousLower':old['lower'],'lower':e['lower'],'upper':e['upper'],'status':e['status'],'subspaceIds':ids})
report={'date':'2026-10-30','catalogueSHA256':hashlib.sha256((ROOT/'data/catalogue.json').read_bytes()).hexdigest(),
        'constructionSHA256':hashlib.sha256((HERE/'bases.json').read_bytes()).hexdigest(),
        'changes':changes,'numberOfImprovements':len(changes),'newRepresentatives':[r['id'] for r in selected],
        'upperBoundsUnchanged':True,'method':'All seven scalar-character Molien families checked at every catalogue qubit number; each exact achieved order propagates to lower orders. Conjugate tetrahedral sectors have equal dimensions. General dihedral scalar-sector dimensions are dominated by the existing first-order family.'}
(HERE/'improvements.json').write_text(json.dumps(report,indent=2)+'\n')
def esc(x):return html.escape(str(x),quote=True)
def math(tex,display=False):return '<span class="math" data-tex="'+esc(tex)+'"'+(' data-display="true"' if display else '')+'></span>'
def poly(coeff):
 x,y=s.symbols('x y');d=len(coeff)-1
 return s.latex(sum(s.sympify(c)*x**r*y**(d-r) for r,c in enumerate(coeff)))
def gf(f):
 num='+'.join('1' if d==0 else 'z^{'+str(d)+'}' for d in f['secondaryDegrees'])
 a,b=f['primaryDegrees']
 return r'\frac{'+num+r'}{(1-z^{'+str(a)+r'})(1-z^{'+str(b)+r'})}'

parts=['''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Molien constructions · Analytical AC subspaces</title><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script defer src="../../vendor/katex/katex.min.js"></script><style>
:root{color-scheme:light;--paper:#f5f4ee;--ink:#27302e;--muted:#66726c;--line:#d6dacf;--purple:#675283;--teal:#245c50}*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:16px/1.8 system-ui,-apple-system,sans-serif}main{max-width:1080px;margin:36px auto;padding:0 28px 50px}h1,h2,h3{font-family:Georgia,serif;font-weight:400;line-height:1.25}h1{font-size:clamp(32px,4vw,46px);letter-spacing:-1px}h2{font-size:29px;margin-top:42px}h3{font-size:22px}a{color:var(--teal);text-underline-offset:4px}.eyebrow{font-size:11px;letter-spacing:.15em;text-transform:uppercase;color:var(--purple);font-weight:700}.lead{font-size:18px;color:#52645a}.result{padding:18px 24px;border:1px solid #c2d3c5;border-left:4px solid var(--teal);background:#eaf1e8;margin:26px 0}.result p:last-child{margin-bottom:0}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;background:#fffef9}th,td{text-align:left;padding:12px 16px;border-bottom:1px solid var(--line);vertical-align:top}th{font-size:12px;color:var(--purple)}.math[data-display=true]{display:block;overflow-x:auto;overflow-y:hidden;max-width:100%;padding:8px 0}.katex-display{margin:.4em 0}.basis{background:#fffef9;border:1px solid var(--line);padding:18px 22px;margin:18px 0}.basis h3{margin:0 0 12px}.note{font-size:14px;color:var(--muted)}details{border-top:1px solid var(--line);padding:13px 0;margin-top:12px}summary{cursor:pointer;color:var(--purple);font-weight:600}.downloads{display:flex;flex-wrap:wrap;gap:12px 24px}footer{border-top:1px solid var(--line);margin-top:40px;padding-top:16px;font-size:13px;color:var(--muted)}:focus-visible{outline:2px solid var(--purple);outline-offset:4px}@media(max-width:600px){main{padding:0 18px}body{font-size:15px}th,td{padding:10px 12px}.basis{padding:15px}.basis .math{font-size:14px}}
</style><script>document.addEventListener('DOMContentLoaded',()=>document.querySelectorAll('[data-tex]').forEach(el=>{if(window.katex)el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==='true',throwOnError:false,trust:false,strict:'ignore'});}));</script></head><body><main><p class="eyebrow">Anticoherent subspaces · Analytical constructions</p><a href="../../index.html">← Return to the catalogue</a><h1>Molien series become explicit AC bases</h1><p class="lead">Invariant and single-character polynomial spaces provide exact subspaces. Every superposition is anticoherent, and all cross conditions are checked exactly.</p>''']
parts.append('<div class="result"><strong>'+str(len(changes))+' stronger rigorous lower bounds</strong><p>'+str(len(source['subspaces']))+' nonzero constructions across seven symmetry families were checked. The catalogue adds '+str(len(selected))+' representatives that improve its existing bounds. Upper bounds are unchanged; these six maxima remain open.</p></div>')
parts.append('<h2>Every improved catalogue entry</h2><div class="table-wrap"><table><thead><tr><th>Qubits</th><th>AC order</th><th>Previous interval</th><th>New interval</th></tr></thead><tbody>')
for e in changes:
 url='../../index.html#n='+str(e['n'])+'&amp;t='+str(e['t'])+'&amp;subspace='+esc(e['subspaceIds'][0])
 parts.append('<tr><td><a href="'+url+'">'+str(e['n'])+'</a></td><td>'+str(e['t'])+'</td><td>'+str(e['previousLower'])+'–'+str(e['upper'])+'</td><td><strong>'+str(e['lower'])+'–'+str(e['upper'])+'</strong></td></tr>')
parts.append('</tbody></table></div><p class="note">A construction at order '+math('t')+' also supplies all orders below '+math('t')+'. Numerical reports are not used to raise these lower bounds.</p>')
parts.append('<h2>Why the entire span is anticoherent</h2><p>Fix a one-dimensional unitary character '+math(r'\chi')+' of a finite spin-rotation group. Define</p>'+math(r'\mathcal C_{G,\chi}(N)=\{\psi:U_N(g)\psi=\chi(g)\psi\quad\forall g\in G\}.',True))
parts.append('<p>For two vectors in this same space, the character phases cancel in each matrix element. Hence, for every spin tensor '+math('T_{Lq}')+',</p>'+math(r'\langle u|T_{Lq}|v\rangle=\left\langle u\left|\frac1{|G|}\sum_{g\in G}U_N(g)^\dagger T_{Lq}U_N(g)\right|v\right\rangle.',True))
parts.append('<p>The averaged tensor vanishes for ranks 1 through 1, 2, 3 and 5 for the quaternion, tetrahedral, octahedral and icosahedral groups respectively. Therefore every low-rank compressed tensor vanishes, including off-diagonal basis pairs. Different characters are kept separate: their direct sum is not automatically AC.</p><p>The sign character of the octahedral group is essential at 22 qubits, where its fixed space is zero but its sign space contains a third-order AC state.</p>')
parts.append('<h2>The seven Molien families</h2><p>The coefficient of '+math('z^N')+' is the dimension at qubit number '+math('N')+'. The linked OEIS sequences index the even coefficients by spin '+math('m=N/2')+'. These dimensions describe the specified symmetry spaces, not unrestricted maxima.</p><div class="table-wrap"><table><thead><tr><th>Family</th><th>Guaranteed order</th><th>Generating function</th><th>Sequence</th></tr></thead><tbody>')
for f in source['families']:
 seq='<a href="https://oeis.org/'+esc(f['oeis'])+'">'+esc(f['oeis'])+'</a>' if f.get('oeis') else 'Character sector'
 parts.append('<tr><td>'+esc(f['name'])+'</td><td>'+str(f['guaranteedOrder'])+'</td><td>'+math(gf(f))+'</td><td>'+seq+'</td></tr>')
parts.append('</tbody></table></div><p class="note">Conjugate tetrahedral characters have equal dimensions; one is exported. The nontrivial quaternion characters are equivalent under a change of axes. Binary dihedral scalar-character spaces have dimension at most '+math(r'\lfloor N/4\rfloor+1')+' and cannot improve the existing first-order construction '+math(r'\lfloor N/3\rfloor+1')+' for '+math(r'N\ge4')+'. Odd qubit numbers have no nonzero scalar-character space in these polyhedral families.</p>')
parts.append('<h2>Polynomial generators and analytical bases</h2><p>Each family is a free module over two primary invariants '+math('A,B')+', with two secondaries '+math('S_0,S_1')+'. Its homogeneous basis consists of the monomials '+math('A^aB^bS_c')+' of total degree '+math('N')+'. The secondaries account for the algebraic relation among the classical generators; counting all unrestricted products would overcount.</p>')
for f in source['families']:
 parts.append('<details id="family-'+esc(f['id'])+'"><summary>'+esc(f['name'])+' · exact generators</summary>')
 for label,p in zip(['A','B','S_0','S_1'],f['primaryPolynomials']+f['secondaryPolynomials']):
  parts.append(math(label+'='+poly(p),True))
 parts.append('</details>')
parts.append('<p>For '+math(r'f(x,y)=\sum_{r=0}^{N}f_r x^ry^{N-r}')+', use the exact weighted inner product</p>'+math(r'\langle f,g\rangle_N=\sum_{r=0}^{N}\frac{\overline{f_r}g_r}{\binom Nr}.',True)+'<p>Exact Gram–Schmidt gives orthogonal polynomials '+math('q_a')+'. Their normalized Dicke states are</p>'+math(r'|u_a\rangle=\frac1{\sqrt{\langle q_a,q_a\rangle_N}}\sum_{r=0}^N\frac{(-1)^r(q_a)_r}{\sqrt{\binom Nr}}|D_N^{(r)}\rangle.',True))
families={f['id']:f for f in source['families']}
for r in selected:
 st=records[r['id']];f=families[st['group']];monomials=[]
 for term in st['sourceMonomials']:
  a,b=term['primaryExponents'];c=term['secondaryIndex']
  monomials.append(('A^{'+str(a)+'}' if a else '')+('B^{'+str(b)+'}' if b else '')+'S_{'+str(c)+'}')
 parts.append('<article class="basis"><h3>'+str(st['n'])+' qubits · dimension '+str(st['dimension'])+' · AC order '+str(st['t'])+'</h3><p>'+esc(f['name'])+'; use its generators above.</p>'+math(r'\mathcal C=\operatorname{span}\{'+',\,'.join(monomials)+r'\}',True)+'<p><a href="../../index.html#n='+str(st['n'])+'&amp;t='+str(st['t'])+'&amp;subspace='+esc(st['id'])+'">Inspect the normalized analytical Dicke basis and Majorana spheres ↗</a></p></article>')
parts.append('<h2>Exact verification and reproducibility</h2><p>The construction checker verifies normalization, orthogonality and every ordered reduced cross matrix. A separate verifier reconstructs the polynomial ranks, rotation-character dimensions and all compressed spin moments, without importing the construction code. Each space is tested above its group-guaranteed order until an exact failure or the catalogue limit; this is a claim about this space, not a higher-order nonexistence proof.</p>'+math(r'R_t(|u_a\rangle\langle u_b|)=\delta_{ab}\frac{I_{\mathrm{Sym}^t}}{t+1}.',True)+'<div class="downloads"><a href="proof.tex">Complete analytical proof ↓</a><a href="bases.json">All polynomial and Dicke bases ↓</a><a href="construct.py">Reproduce the constructions ↓</a><a href="verification.json">Construction verification ↓</a><a href="verify_independently.py">Independent verifier ↓</a><a href="independent-verification.json">Independent verification ↓</a><a href="improvements.json">All lower-bound changes ↓</a></div>')
parts.append('<h2>Sources and scope</h2><p>The construction uses the invariant-polynomial method of <a href="https://arxiv.org/html/1610.05841v2#S3">Pereira and Paul-Paddock, Anticoherent Subspaces (2017), Section 3</a>, and the Molien series linked above. The catalogue supplies the explicit orthonormal bases and exact finite verifications for its listed qubit numbers. No claim of a new invariant-theory theorem or of global maximality is made.</p><footer>30 October 2026 · <a href="../../index.html">Return to the AC-subspace catalogue ↗</a></footer></main></body></html>')
(HERE/'index.html').write_text('\n'.join(parts)+'\n')
print(json.dumps({'improvements':len(changes),'representatives':len(selected),'changes':changes},indent=2))
