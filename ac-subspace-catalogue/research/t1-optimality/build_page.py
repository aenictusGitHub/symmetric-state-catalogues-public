#!/usr/bin/env python3
"""Build the first-order optimality proof page and its historical change report."""
from pathlib import Path
import json,html,re,hashlib
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def tex(s,display=False):
 return '<div class="math" data-display="true" data-tex="'+html.escape(s,quote=True)+'"></div>' if display else '<span class="math" data-tex="'+html.escape(s,quote=True)+'"></span>'
def main():
 data=json.loads((ROOT/'data/catalogue.json').read_text())
 old=json.loads((HERE/'baseline-entries.json').read_text())
 prior={e['id']:e for e in old}
 changes=[]
 for e in data['entries']:
  b=prior[e['id']]
  if (e['lower'],e['upper'])!=(b['lower'],b['upper']):
   assert e['t']==1 and e['lower']==e['upper']==e['n']//3+1
   changes.append({'n':e['n'],'t':1,'previousLower':b['lower'],'previousUpper':b['upper'],'lower':e['lower'],'upper':e['upper'],'status':e['status']})
 assert len(changes)==11
 report={'theorem':'K(N,1)=floor(N/3)+1 for all N>=4','proof':'Auxiliary spin-1 inertia and matching analytical construction','improvements':changes,'catalogueSha256':hashlib.sha256((ROOT/'data/catalogue.json').read_bytes()).hexdigest()}
 (HERE/'improvements.json').write_text(json.dumps(report,indent=2)+'\n')
 css=re.search(r'<style>(.*?)</style>',(ROOT/'research/ancilla-inertia/index.html').read_text(),re.S).group(1)
 rows=''.join('<tr><td><a href="../../index.html#n='+str(c['n'])+'&amp;t=1">'+str(c['n'])+'</a></td><td>'+str(c['previousLower'])+'–'+str(c['previousUpper'])+'</td><td>'+str(c['lower'])+' · proved maximum</td></tr>' for c in changes)
 body='''<p class="eyebrow">Anticoherent subspaces · General optimality</p>
<a class="back-link" href="../../index.html#n=14&amp;t=1">← Return to the catalogue</a>
<h1>The first-order maximum, for every qubit number</h1>
<p class="lead">An auxiliary spin-1 system supplies the missing upper bound. The analytical reflection and central blocks attain it, settling the first-order sequence for all sizes.</p>
<div class="result"><strong>General theorem</strong>'''+tex(r'\boxed{K(N,1)=\left\lfloor\frac N3\right\rfloor+1\qquad(N\ge4)}',True)+'''<p>Here '''+tex(r'K(N,1)')+''' is the largest dimension of a complex first-order anticoherent subspace of '''+tex(r'\operatorname{Sym}^{N}(\mathbb C^2)')+'''. Every normalized vector in the subspace must be anticoherent.</p></div>
<h2>1. Zero compression survives the auxiliary spin</h2>
<p>Write '''+tex(r'j=N/2')+''' and let '''+tex(r'P')+''' project onto a first-order AC subspace '''+tex(r'\mathcal C')+''' of dimension '''+tex(r'k')+'''. Polarization of the vanishing spin expectations gives the full matrix identities</p>'''+tex(r'PJ_xP=PJ_yP=PJ_zP=0.',True)+'''
<p>Let '''+tex(r'S_x,S_y,S_z')+''' be spin-1 generators on a three-dimensional auxiliary space and define</p>'''+tex(r'W=J_x\otimes S_x+J_y\otimes S_y+J_z\otimes S_z,\qquad\Pi=P\otimes I_3.',True)+'''
<p>Then</p>'''+tex(r'\Pi W\Pi=\sum_{a=x,y,z}(PJ_aP)\otimes S_a=0.',True)+'''
<p>Thus the Hermitian form defined by '''+tex(r'W')+''' vanishes on the whole '''+tex(r'3k')+'''-dimensional space '''+tex(r'\mathcal C\otimes\mathbb C^3')+''', including entangled vectors in that tensor product.</p>
<h2>2. Three eigenvalues, exactly</h2>
<p>For '''+tex(r'N\ge2')+''', angular-momentum addition gives</p>'''+tex(r'V_j\otimes V_1=V_{j+1}\oplus V_j\oplus V_{j-1},\qquad W=\frac12\bigl((\mathbf J+\mathbf S)^2-j(j+1)I-2I\bigr).',True)+'''
<div class="table-wrap"><table class="classification"><thead><tr><th>Total spin</th><th>Eigenvalue of '''+tex('W')+'''</th><th>Multiplicity</th></tr></thead><tbody>
<tr><td>'''+tex(r'j+1')+'''</td><td>'''+tex('j>0')+'''</td><td>'''+tex('N+3')+'''</td></tr>
<tr><td>'''+tex('j')+'''</td><td>'''+tex('-1')+'''</td><td>'''+tex('N+1')+'''</td></tr>
<tr><td>'''+tex(r'j-1')+'''</td><td>'''+tex(r'-(j+1)')+'''</td><td>'''+tex('N-1')+'''</td></tr>
</tbody></table></div><p>There are exactly '''+tex(r'N+3')+''' positive directions and '''+tex(r'2N')+''' negative directions, with no zero eigenvalues.</p>
<h2>3. Projection proves the upper bound</h2>
<p>Let '''+tex(r'Q_+')+''' be the spectral projection onto the positive eigenspace. Its restriction to '''+tex(r'\mathcal C\otimes\mathbb C^3')+''' is injective. Indeed, a nonzero vector '''+tex(r'v')+''' in its kernel would lie entirely in the negative eigenspace, so '''+tex(r'\langle v,Wv\rangle<0')+'''. Zero compression instead gives '''+tex(r'\langle v,Wv\rangle=0')+''', a contradiction.</p>'''+tex(r'3k=\dim(\mathcal C\otimes\mathbb C^3)\le\dim\operatorname{ran}Q_+=N+3.',True)+'''
<p>Consequently</p>'''+tex(r'K(N,1)\le\left\lfloor\frac{N+3}{3}\right\rfloor=\left\lfloor\frac N3\right\rfloor+1.',True)+'''
<p>Projection into the negative eigenspace also gives '''+tex(r'3k\le2N')+'''. The argument allows arbitrary complex subspaces and imposes no Dicke support, reality, or constellation symmetry.</p>
<h2>4. The analytical construction attains the bound</h2>
<p>The <a href="../analytical-constructions/index.html#blocks">explicit reflection-pair and central blocks</a> give orthonormal bases of dimension '''+tex(r'\lfloor N/3\rfloor+1')+''' for every '''+tex(r'N\ge4')+'''. Their full compressed spin matrices vanish. The support gaps remove inter-block terms; exact radical phases cancel the remaining central-block terms. This is a general construction covering all six residue classes of '''+tex(r'N\bmod6')+'''.</p>
<p>The matching lower and upper bounds establish equality for every size, beyond the finite range shown in the catalogue.</p>
<h2>Small sizes and the integer sequence</h2>'''+tex(r'K(1,1)=0,\qquad K(2,1)=K(3,1)=1.',True)+'''
<p>The one-qubit cross reductions supply '''+tex(r'2k')+''' orthonormal conditional vectors in a space of dimension '''+tex('N')+''', so '''+tex(r'2k\le N')+'''. This settles these small sizes with the two-qubit middle Dicke state and the three-qubit GHZ state as witnesses.</p>
<p>Using the zero-based convention '''+tex(r'a(N)=\lfloor N/3\rfloor+1')+''', the maximum agrees with <a href="https://oeis.org/A008620">OEIS A008620</a> for every '''+tex(r'N\ge4')+'''. The sequence identification follows from the theorem; an OEIS match alone would not prove maximality.</p>
<h2>Eleven newly settled catalogue entries</h2>
<div class="table-wrap"><table class="classification"><thead><tr><th>Qubits</th><th>Previous interval</th><th>Exact maximum at order 1</th></tr></thead><tbody>'''+rows+'''</tbody></table></div>
<p>All first-order cells are now settled. Higher-order values and all analytical representatives are unchanged. The earlier finite-size exterior-power certificates remain available as independent proofs.</p>
<h2>Proof and reproducibility</h2>
<p><a href="proof.tex">Complete LaTeX proof ↓</a> · <a href="verify.py">Exact spectral verifier ↓</a> · <a href="verification.json">Spectral verification ↓</a> · <a href="verify_independently.py">Independent physical-matrix verifier ↓</a> · <a href="independent-verification.json">Independent verification ↓</a> · <a href="improvements.json">All changed bounds ↓</a></p>
<p>The angular-momentum decomposition and the injective-projection argument prove the statement for all sizes. Finite exact matrix checks validate conventions and implementation. No claim of literature priority is made.</p>
<footer><a href="../analytical-constructions/index.html">Optimal analytical bases ↗</a></footer>'''
 head='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>First-order AC subspaces: the exact maximum for every N</title><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script defer src="../../vendor/katex/katex.min.js"></script><style>'+css+'</style><script>document.addEventListener("DOMContentLoaded",()=>document.querySelectorAll("[data-tex]").forEach(el=>{if(window.katex)el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==="true",throwOnError:false,strict:"ignore",trust:false});}));</script></head><body><main>'
 (HERE/'index.html').write_text(head+body+'</main></body></html>\n')
 print(json.dumps({'improvedEntries':len(changes),'validation':data['validation']},indent=2))
if __name__=='__main__':main()
