/* Static, dependency-free catalogue UI; KaTeX is vendored for local typesetting. */
(() => {
  'use strict';
  const data = window.SECTOR_CATALOGUE;
  if (!data || !Array.isArray(data.entries)) return;
  const $ = id => document.getElementById(id);
  const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const labels = {proven:'Proved / analytic',reported:'Conjectural optimum',numerical:'Numerical candidate',bound:'Bound only',unrecorded:'Not catalogued'};
  const codes = {proven:'P',reported:'C',numerical:'N',bound:'B',unrecorded:'—'};
  const math = (tex, display=false) => {
    try { return window.katex.renderToString(tex, {displayMode:display,throwOnError:false,trust:false,output:'htmlAndMathml'}); }
    catch { return `<span class="math">${escape(tex)}</span>`; }
  };
  const rich = value => String(value ?? '').split(/(\$[^$]+\$)/g).map(part => part.startsWith('$') && part.endsWith('$') ? math(part.slice(1,-1)) : escape(part)).join('');
  const sourceById = new Map(data.sources.map(s => [s.id,s]));
  const recordByKey = new Map(data.entries.map(e => [`${e.n}-${e.k}-${e.mode}`,e]));
  let mode = 'min', selectedN = 6, selectedK = 3, evidence = 'all';
  const numbers = data.qubitNumbers || [...new Set(data.entries.map(e => e.n))].sort((a,b) => a-b);
  const maxK = Math.max(...numbers);
  const selected = () => recordByKey.get(`${selectedN}-${selectedK}-${mode}`);
  const sourceRef = (ref) => {
    const source = sourceById.get(ref.sourceId);
    if (!source) return '';
    const location = [ref.page ? `PDF p. ${ref.page}` : '',ref.location].filter(Boolean).join(' · ');
    const url=ref.url || `${source.pdf || source.href || source.url}${ref.page ? '#page='+ref.page : ref.anchor ? '#'+ref.anchor : ''}`;
    return `<a class="source-link" href="${escape(url)}" target="_blank" rel="noopener">${escape(source.shortTitle || source.title)}${location ? ` — ${escape(location)}` : ''} ↗</a>`;
  };
  function readHash() {
    const params = new URLSearchParams(location.hash.slice(1));
    const n = Number(params.get('n')), k = Number(params.get('k'));
    if (numbers.includes(n) && Number.isInteger(k) && k>=1 && k<=n) {selectedN=n;selectedK=k;}
    if (['min','max'].includes(params.get('mode'))) mode=params.get('mode');
  }
  function updateHash() { history.replaceState(null,'',`#n=${selectedN}&k=${selectedK}&mode=${mode}`); }
  function shortValue(entry) {
    if (!entry) return '—';
    if (entry.cellLabel) return entry.cellLabel;
    if (String(entry.value).length <= 8) return entry.value;
    return Number.isFinite(entry.decimal) ? Number(entry.decimal.toPrecision(5)).toString() : 'exact';
  }
  function renderTable() {
    $('mode-min').setAttribute('aria-pressed',String(mode==='min'));
    $('mode-max').setAttribute('aria-pressed',String(mode==='max'));
    let html = `<table class="sector-grid"><caption class="visually-hidden">${mode==='min'?'Minimum':'Maximum'} reduced purities for pure symmetric states. Rows: qubits n; columns: reduction size t. P: proved / analytic; C: conjectural optimum; N: numerical candidate; B: bound only.</caption><thead><tr><th scope="col">${math('n\\backslash t')}</th>`;
    for(let k=1;k<=maxK;k++) html += `<th scope="col">${k}</th>`;
    html+='</tr></thead><tbody>';
    for (const n of numbers) {
      html+=`<tr><th scope="row">${n}</th>`;
      for(let k=1;k<=maxK;k++) {
        if(k>n) { html+='<td class="not-applicable" aria-label="Reduction exceeds number of qubits"></td>'; continue; }
        const entry=recordByKey.get(`${n}-${k}-${mode}`), status=entry?.status || 'unrecorded';
        const isSelected=selectedN===n&&selectedK===k, muted=evidence!=='all'&&status!==evidence;
        const accessible=`${n} qubits, reduction ${k}, ${mode==='min'?'minimum':'maximum'}: ${entry?entry.value:'not catalogued'}, ${labels[status]}`;
        html+=`<td><button type="button" class="sector-cell ${status}${isSelected?' selected':''}${muted?' filtered-out':''}" data-n="${n}" data-k="${k}" aria-pressed="${isSelected}" aria-label="${escape(accessible)}" title="${escape(accessible)}"><span>${escape(shortValue(entry))}</span><small aria-hidden="true">${codes[status]}</small></button></td>`;
      }
      html+='</tr>';
    }
    $('sector-table').innerHTML=html+'</tbody></table>';
    const entries=data.entries.filter(e=>e.mode===mode&&(evidence==='all'||e.status===evidence));
    $('result-count').textContent=`${entries.length} ${evidence==='all'?'catalogued':'matching'} results · n = ${numbers[0]}–${numbers.at(-1)}`;
  }
  function renderDetail() {
    window.MajoranaViewer?.close();
    const entry=selected();
    const majorana=entry && window.MAJORANA_DATA?.records?.[entry.id];
    let html=`<div class="detail-top"><p class="eyebrow">${entry && ['reported','numerical'].includes(entry.status)?'Candidate ':''}${mode==='min'?'minimum':'maximum'} / selected result</p><a href="#catalogue" class="detail-back">↑ Table</a></div><h3>${selectedN} qubits · ${math(`P_{${selectedK}}`)}</h3>`;
    if (!entry) {
      $('result-detail').innerHTML=html+'<span class="badge unrecorded">Not catalogued</span><p class="detail-note">No extremum for this cell is recorded in the supplied source tables or the general results used here. This is a gap in this catalogue, not a claim that the problem is open.</p>';
      return;
    }
    html+=`<span class="badge ${entry.status}">${labels[entry.status]}</span>`;
    html+=`<div class="detail-value${entry.value.length>16?' long-value':''}">${entry.valueLatex?math(entry.valueLatex,true):escape(entry.value)}</div>`;
    if(entry.decimal!==null && entry.decimal!==undefined && String(entry.value)!==String(entry.decimal)) html+=`<p class="decimal-value">≈ ${Number(entry.decimal.toPrecision(12))}</p>`;
    if(entry.exactValue) html+=`<div class="bound-callout"><h4>Exact attained value</h4><p>${math(entry.exactValue.formulaLatex)}</p><p>${escape(entry.exactValue.summary)}</p><a href="${escape(entry.exactValue.url)}">Exact expression & proof ↗</a></div>`;
    if(entry.certifiedBound && entry.status!=='proven') html+=`<div class="bound-callout"><h4>Certified ${mode==='min'?'lower':'upper'} bound</h4><p>${math(`P_{${entry.k}} ${mode==='min'?'\\ge':'\\le'} ${entry.certifiedBound.valueLatex}`)}</p><a href="${escape(entry.certifiedBound.url || entry.proofUrl || '#definitions')}">View the remaining gap & certificate ↗</a></div>`;
    html+=`<p class="evidence-description">${escape(entry.evidence || (entry.status==='proven'?'Global optimality is established in the cited source.':entry.status==='reported'?'An exact candidate is reported; a global proof is not supplied for this entry.':entry.status==='numerical'?'Numerical optimization supplies this candidate; global optimality is not proved.':'The source gives a bound; attainability is not established.'))}</p>`;
    if(entry.proofUrl) html+=`<p><a href="${escape(entry.proofUrl)}">Analytical bound and verification ↗</a></p>`;
    if(entry.state || entry.stateLatex) {
      html+='<h4>Representative state</h4>';
      if(entry.stateName) html+=`<p class="state-name">${escape(entry.stateName)}</p>`;
      html+=window.MajoranaViewer?.groupSummary?.(majorana)||'';
      html+=`<div class="state-vector">${entry.stateLatex?math(entry.stateLatex,true):rich(entry.state)}</div>`;
    }
    if(entry.parameters) html+=`<p>${rich(entry.parameters)}</p>`;
    if(majorana?.points?.length)html+='<button id="view-majorana" class="view-majorana-button" type="button" aria-haspopup="dialog" aria-controls="majorana-viewer">View Majorana sphere ↗</button>';
    if(entry.algebraicUrl) html+=`<p><a href="${escape(entry.algebraicUrl)}">Exact algebraic coefficients and derivation ↗</a></p>`;
    if(entry.coefficients) html+=`<details class="coefficient-detail"><summary>Dicke coefficients</summary><p class="muted">Nonzero amplitudes; full precision is in the downloadable record.</p><table><thead><tr><th>r</th><th>cᵣ</th></tr></thead><tbody>${entry.coefficients.map((z,r)=>({z:Array.isArray(z)?z:[z,0],r})).filter(({z})=>Math.hypot(...z)>1e-12).map(({z,r})=>`<tr><td>${r}</td><td>${Number(z[0].toPrecision(9))}${Math.abs(z[1])>1e-12?`${z[1]>=0?'+':''}${Number(z[1].toPrecision(9))}i`:''}</td></tr>`).join('')}</tbody></table></details>`;
    if(entry.notes) html+=entry.notes.length>240?`<details class="proof-notes"><summary>Proof and scope notes</summary><p class="detail-note">${rich(entry.notes)}</p></details>`:`<p class="detail-note">${rich(entry.notes)}</p>`;
    if(entry.sectorVector) html+=`<details class="sector-vector"><summary>Full sector vector of this representative</summary><p class="muted">Computed from its Dicke coefficients; entries are ${math('S_1,\\ldots,S_n')}.</p><p class="vector-values">(${entry.sectorVector.map(v=>Number(v.toPrecision(7))).join(', ')})</p></details>`;
    html+='<h4>Source & evidence</h4><ul class="detail-sources">';
    for(const ref of entry.references || []) html+=`<li>${sourceRef(ref)}</li>`;
    html+='</ul>';
    if(entry.certificateUrl) html+=`<p><a class="source-link" href="${escape(entry.certificateUrl)}" target="_blank" rel="noopener">Certificate repository cited by the draft ↗</a></p>`;
    if(entry.status==='proven' && !entry.proofUrl && entry.references?.some(r=>sourceById.get(r.sourceId)?.draft)) html+='<p class="draft-note">Includes a result proved in the supplied working draft. See the source record for its publication status.</p>';
    html+='<div class="detail-actions"><button type="button" id="download-result">Download record ↓</button><button type="button" id="copy-link">Copy link ↗</button></div><p id="action-feedback" role="status"></p>';
    $('result-detail').innerHTML=html;
    $('view-majorana')?.addEventListener('click',()=>window.MajoranaViewer?.open(majorana));
    $('download-result').addEventListener('click',()=>download(JSON.stringify({...entry,majoranaRepresentation:majorana||null},null,2),`symmetric-n${selectedN}-p${selectedK}-${mode}.json`,'application/json'));
    $('copy-link').addEventListener('click',async()=>{
      updateHash();
      try {await Promise.race([navigator.clipboard.writeText(location.href),new Promise((_,reject)=>setTimeout(()=>reject(new Error('Clipboard unavailable')),1500))]);$('action-feedback').textContent='Link copied.';}
      catch {$('action-feedback').innerHTML=`Share this address: <input aria-label="Result link" readonly value="${escape(location.href)}">`;}
    });
  }
  function download(content,filename,type) {
    const url=URL.createObjectURL(new Blob([content],{type}));
    const a=document.createElement('a');a.href=url;a.download=filename;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  }
  function renderReferences() {
    $('coverage-note').innerHTML=rich(data.coverageNote);
    $('constructions-list').innerHTML=(data.constructions||[]).map(c=>`<article><h3>${escape(c.title)}</h3><p>${rich(c.text)}</p>${c.formula?`<div class="construction-formula">${math(c.formula,true)}</div>`:''}<p class="construction-source">${(c.references||[]).map(sourceRef).join('<br>')}</p></article>`).join('');
    $('definitions-content').innerHTML=(data.definitions||[]).map(d=>`<article><h3>${escape(d.title)}</h3><p>${rich(d.text)}</p>${d.formula?`<div class="construction-formula">${math(d.formula,true)}</div>`:''}${d.references?`<p>${d.references.map(sourceRef).join('<br>')}</p>`:''}</article>`).join('');
    $('source-list').innerHTML=`<p class="source-intro">Local computation and analytical witness verification. Numerical candidates remain explicitly distinguished from proofs.</p>`+data.sources.map(s=>`<article id="source-${escape(s.id)}"><div class="source-meta"><span>${escape(s.role)}</span><span>${escape(s.year)}${s.draft?' · Working manuscript':''}</span></div><h3>${escape(s.title)}</h3><p>${escape(s.authors)}</p><p>${rich(s.summary)}</p>${s.notes?`<p class="source-notes">${rich(s.notes)}</p>`:''}<p class="source-actions">${s.pdf?`<a href="${escape(s.pdf)}" target="_blank" rel="noopener">Source PDF ↗</a>`:''}${s.href?`<a href="${escape(s.href)}">Proofs & computations ↗</a>`:''}${s.doi?`<a href="https://doi.org/${escape(s.doi)}" target="_blank" rel="noopener">DOI ↗</a>`:''}${s.url?`<a href="${escape(s.url)}" target="_blank" rel="noopener">Public source ↗</a>`:''}<span>${escape(s.filename)}</span></p></article>`).join('');
  }
  $('sector-table').addEventListener('click',event=>{
    const button=event.target.closest('button[data-n]');if(!button)return;
    selectedN=Number(button.dataset.n);selectedK=Number(button.dataset.k);updateHash();
    document.querySelectorAll('.sector-cell.selected').forEach(el=>{el.classList.remove('selected');el.setAttribute('aria-pressed','false');});
    button.classList.add('selected');button.setAttribute('aria-pressed','true');renderDetail();
    if(matchMedia('(max-width:800px)').matches) $('result-detail').scrollIntoView({behavior:matchMedia('(prefers-reduced-motion:reduce)').matches?'instant':'smooth',block:'start'});
  });
  $('sector-table').addEventListener('keydown',event=>{
    const button=event.target.closest('button[data-n]');if(!button || !['ArrowUp','ArrowDown','ArrowLeft','ArrowRight'].includes(event.key))return;
    event.preventDefault();let n=Number(button.dataset.n),k=Number(button.dataset.k);
    if(event.key==='ArrowUp')n--;if(event.key==='ArrowDown')n++;if(event.key==='ArrowLeft')k--;if(event.key==='ArrowRight')k++;
    if(numbers.includes(n)&&k>=1&&k<=n) document.querySelector(`button[data-n="${n}"][data-k="${k}"]`)?.focus();
  });
  ['min','max'].forEach(value=>$(`mode-${value}`).addEventListener('click',()=>{mode=value;updateHash();renderTable();renderDetail();}));
  $('status-filter').addEventListener('change',event=>{evidence=event.target.value;renderTable();});
  window.addEventListener('hashchange',()=>{readHash();renderTable();renderDetail();});
  readHash();renderTable();renderDetail();renderReferences();
})();
