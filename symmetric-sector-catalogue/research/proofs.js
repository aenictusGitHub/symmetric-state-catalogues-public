(() => {
  'use strict';
  const data=window.SECTOR_CATALOGUE;
  const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const math=(s,display=false)=>katex.renderToString(s,{displayMode:display,throwOnError:false,trust:false});
  const frac=s=>{const [a,b]=s.split('/');return b?`\\frac{${a}}{${b}}`:a;};
  function term(name) {
    if(/^p_\d+$/.test(name))return `p_{${name.slice(2)}}`;
    if(/^a_\d+$/.test(name))return `a_{${name.slice(2)}}`;
    if(name.startsWith('even_reduced_upper_'))return `E_{${name.split('_').at(-1)}}`;
    const purity=name.match(/^r_(\d+)-1\/(\d+)$/);
    if(purity)return `\\left(r_{${purity[1]}}-\\frac1{${purity[2]}}\\right)`;
    return `\\mathrm{g}_{\\text{${name.replaceAll('_',' ').replace(/[{}\\]/g,'')}}}`;
  }
  function certificate(e,c) {
    const bound=e.status==='proven'?e.valueLatex:e.certifiedBound.valueLatex;
    const lhs=e.mode==='min'?`S_{${e.k}}-${bound}`:`${bound}-S_{${e.k}}`;
    const terms=Object.entries(c).map(([key,value])=>(value==='1'?'':frac(value)+'\\,')+term(key));
    return math(`\\begin{aligned}${lhs}&=${terms.slice(0,2).join('+')}${terms.length>2?'\\\\&\\quad+'+terms.slice(2,4).join('+'):''}${terms.length>4?'\\\\&\\quad+'+terms.slice(4,6).join('+'):''}${terms.length>6?'\\\\&\\quad+'+terms.slice(6).join('+'):''}\\ge0.\\end{aligned}`,true);
  }
  document.querySelectorAll('[data-tex]').forEach(el=>el.innerHTML=math(el.dataset.tex));
  document.querySelectorAll('[data-display-tex]').forEach(el=>el.innerHTML=math(el.dataset.displayTex,true));
  const sharp=data.entries.filter(e=>e.proofCertificate);
  document.getElementById('sharp-results').innerHTML=sharp.map(e=>`<article class="proof-card" id="${e.id}"><div class="proof-card-heading"><h3>${e.n} qubits · ${e.mode==='min'?'minimum':'maximum'} ${math(`S_{${e.k}}`)} = ${math(e.valueLatex)}</h3><span class="badge proven">Exact proof</span></div><div class="certificate-equation">${certificate(e,e.proofCertificate)}</div><p>Attained by</p><div class="witness-equation">${math(e.stateLatex,true)}</div>${e.parameters?`<p>${e.parameters.split(/(\$[^$]+\$)/g).map(s=>s.startsWith('$')?math(s.slice(1,-1)):esc(s)).join('')}</p>`:''}<p class="muted">${esc(e.n>=10&&e.k===4&&e.mode==='min'?'The source construction cancels the even multipoles through rank four; normalization and equality were independently verified symbolically.':'The identity is checked with exact rational arithmetic and the displayed state with exact symbolic arithmetic.')}</p><a href="../index.html#n=${e.n}&k=${e.k}&mode=${e.mode}">Open catalogue record & coefficients ↗</a></article>`).join('');
  const gaps=data.entries.filter(e=>e.status!=='proven');
  document.querySelector('.gap-table-wrap').hidden=gaps.length===0;
  document.getElementById('gap-summary').textContent=gaps.length ? `${gaps.length} candidate optima remain. Each attained candidate has a separate certified global bound.` : `All ${data.entries.length} catalogue extrema are proved. Historical relaxation bounds and numerical discovery records remain available through their proof pages.`;
  document.getElementById('gap-rows').innerHTML=gaps.map(e=>`<tr><td>${e.n}</td><td>${math(`S_{${e.k}}`)}</td><td>${e.mode}</td><td><a href="#${e.id}">${math(e.valueLatex||e.value)}</a></td><td>${math((e.mode==='min'?'\\ge':'\\le')+e.certifiedBound.valueLatex)}</td><td>${e.status==='numerical'?'Numerical state':'Exact state; optimum conjectural'}</td></tr>`).join('');
  document.getElementById('gap-certificates').innerHTML=gaps.map(e=>`<details class="gap-certificate" id="${e.id}"><summary>${e.n} qubits · ${e.mode} S${e.k} — candidate ${esc(e.value)}</summary><div class="certificate-equation">${e.certifiedBound.certificate?certificate(e,e.certifiedBound.certificate):math(`S_{${e.k}}${e.mode==='min'?'\\ge':'\\le'}${e.certifiedBound.valueLatex}`,true)}</div>${e.certifiedBound.certificate?`<p>All named quantities are nonnegative under the stated assumptions. For terms beyond ${math('p_J')} and ${math('a_L')}, see <a href="PROOFS.md">the full proof assumptions</a> and <a href="certificates/certificate_basis.py">the exact definitions</a>.</p>`:`<p>${esc(e.certifiedBound.method)}. <a href="../${esc(e.certifiedBound.proofUrl)}">Read the exact bound proof and verification ↗</a></p>`}<div class="witness-equation">${math(e.stateLatex,true)}</div><p>${esc(e.notes)}</p>${e.algebraicUrl?`<p><a href="../${esc(e.algebraicUrl)}">Exact algebraic coefficients and derivation ↗</a></p>`:''}<a href="../index.html#n=${e.n}&k=${e.k}&mode=${e.mode}">Open candidate and download coefficients ↗</a></details>`).join('');
  function reveal(){const id=location.hash.slice(1);const target=document.getElementById(id);if(target?.tagName==='DETAILS')target.open=true;}
  reveal();window.addEventListener('hashchange',reveal);
})();
