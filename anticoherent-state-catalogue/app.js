/* Offline catalogue viewer. Scientific records live in data/catalogue.js. */
(() => {
  'use strict';
  const data = window.AC_CATALOGUE;
  const $ = id => document.getElementById(id);
  const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const math = (latex, display=false, cls='') => `<span class="math ${cls}" data-tex="${escape(latex)}" data-display="${display}">${escape(latex)}</span>`;
  const labels = {exists:'Explicit witness', impossible:'Proven impossible', reported_absent:'Reported absent', unknown:'Unresolved'};
  const shortLabels = {exists:'exists', impossible:'impossible', reported_absent:'reported no', unknown:'unresolved'};
  const symbols = {exists:'●', impossible:'×', reported_absent:'?', unknown:'—'};
  if (!data || !Array.isArray(data.entries)) {
    $('state-table').innerHTML = '<div class="empty-state"><span class="empty-mark" aria-hidden="true">◎</span><h3>Source records are not loaded</h3><p>The catalogue data is needed to display verified results.</p></div>';
    $('result-count').textContent = 'Source data pending';
    return;
  }
  const states = new Map((data.states || []).map(s => [s.id,s]));
  const constellation = state => state && window.MAJORANA_DATA?.records?.[state.id];
  const symmetrySummary = (state,compact=false) => window.MajoranaViewer?.groupSummary?.(constellation(state),compact)||'';
  const symmetryLabel = state => {const label=window.MajoranaViewer?.groupLabel?.(constellation(state));return label?` · ${label}`:'';};
  const entries = new Map(data.entries.map(e => [`${e.n}:${e.t}`,e]));
  const proofs = new Map((data.proofs || []).map(p => [p.id,p]));
  let selected = {n:12,t:5,stateId:null};
  let filter = 'all';
  let extended = false;
  let toastTimer;

  function typeset(root=document) {
    root.querySelectorAll('[data-tex]').forEach(el => {
      if (window.katex) el.innerHTML = window.katex.renderToString(el.dataset.tex, {displayMode:el.dataset.display==='true',throwOnError:false,trust:false,strict:'ignore'});
    });
  }
  function numbers() { return [...new Set([...(data.coreNumbers || []), ...(extended ? data.extendedNumbers || [] : [])])].sort((a,b)=>a-b); }
  function orders() { return (data.orders || [1,2,3,4,5,6,7]).filter(t=>extended || t<=6); }
  function getEntry(n,t) { return entries.get(`${n}:${t}`) || {id:`n${n}-t${t}`,n,t,status:'unknown',stateIds:[],notes:'No existence result or exclusion is recorded for this pair.'}; }
  function matches(entry) { return filter==='all' || entry.status===filter; }
  function status(entry) { return Object.hasOwn(labels,entry.status) ? entry.status : 'unknown'; }
  function badge(entry) { const s=status(entry); return `<span class="badge ${s}">${labels[s]}</span>`; }
  function sourceUrl(ref) {
    if (ref.file || ref.filename) return `references/${encodeURIComponent(ref.file || ref.filename)}${ref.page ? '#page='+Number(ref.page) : ''}`;
    if (ref.pdf) return `${ref.pdf.includes('/') ? ref.pdf : 'references/'+encodeURIComponent(ref.pdf)}${ref.page ? '#page='+Number(ref.page) : ''}`;
    return ref.url || '#sources';
  }
  function safeUrl(url) { return /^(?:https?:|\.\.?\/|references\/|#)/i.test(url) ? url : '#sources'; }
  function referenceLinks(refs=[]) {
    const unique = new Map(refs.map(r=>[`${r.file||r.filename||r.url}:${r.page||''}:${r.location||''}`,r]));
    return [...unique.values()].map(ref=>`<li><a class="source-link" href="${escape(safeUrl(sourceUrl(ref)))}"${ref.url && !ref.file ? ' target="_blank" rel="noopener"' : ''}>${escape(ref.title || ref.file || ref.filename || 'Source')}${ref.page ? ` · p. ${escape(ref.page)}` : ''}<span aria-hidden="true"> ↗</span></a>${ref.location ? `<span class="source-location">${escape(ref.location)}</span>` : ''}</li>`).join('');
  }
  function cleanBody(html='') {
    const template = document.createElement('template'); template.innerHTML=html;
    const allowed=new Set(['P','UL','OL','LI','EM','STRONG','CODE','PRE','SPAN','DIV','A','SUB','SUP','BR','DETAILS','SUMMARY','BLOCKQUOTE']);
    [...template.content.querySelectorAll('*')].forEach(el=>{
      if (!allowed.has(el.tagName)) { el.replaceWith(document.createTextNode(el.textContent)); return; }
      [...el.attributes].forEach(attr=>{
        if (!['class','data-tex','data-display','href'].includes(attr.name)) el.removeAttribute(attr.name);
        if (attr.name==='href') el.setAttribute('href',safeUrl(attr.value));
      });
    });
    return template.innerHTML;
  }
  function notify(message) { $('toast').textContent=message; $('toast').classList.add('visible'); clearTimeout(toastTimer); toastTimer=setTimeout(()=>$('toast').classList.remove('visible'),2600); }
  function hashFor() { const p=new URLSearchParams({n:selected.n,t:selected.t}); if(selected.stateId)p.set('state',selected.stateId); return '#'+p.toString(); }
  function updateHash() { const hash=hashFor(); if(location.hash!==hash) history.replaceState(null,'',hash); }
  function readHash() {
    const params=new URLSearchParams(location.hash.slice(1)); const n=Number(params.get('n')),t=Number(params.get('t'));
    if ([...(data.coreNumbers||[]),...(data.extendedNumbers||[])].includes(n) && (data.orders||[]).includes(t) && entries.has(`${n}:${t}`)) {
      selected={n,t,stateId:params.get('state')};
      if ((data.extendedNumbers||[]).includes(n) || t>6) extended=true;
    }
  }
  function select(n,t,stateId=null,{focus=false}={}) {
    selected={n:Number(n),t:Number(t),stateId};
    if ((data.extendedNumbers||[]).includes(selected.n) || selected.t>6) extended=true;
    if (!matches(getEntry(selected.n,selected.t))) filter='all';
    syncControls(); renderMatrix(); renderDetail(); updateHash();
    if(focus) $('state-table').querySelector(`[data-n="${selected.n}"][data-t="${selected.t}"]`)?.focus({preventScroll:true});
  }
  function syncControls() { $('status-filter').value=filter; $('extended-toggle').checked=extended; }
  function ensureSelection() {
    const ns=numbers(),ts=orders();
    if (!ns.includes(selected.n) || !ts.includes(selected.t) || !matches(getEntry(selected.n,selected.t))) {
      const first=ns.flatMap(n=>ts.map(t=>entries.get(`${n}:${t}`)).filter(Boolean)).find(matches);
      if(first)selected={n:first.n,t:first.t,stateId:null};
    }
  }
  function renderMatrix() {
    const ns=numbers(),ts=orders(); let matching=0;
    const head=`<thead><tr><th class="corner" scope="col">${math('n\\backslash t')}</th>${ts.map(t=>`<th scope="col">${t}</th>`).join('')}</tr></thead>`;
    const body=ns.map(n=>`<tr><th scope="row">${n}</th>${ts.map(t=>{
      if(t>n || !entries.has(`${n}:${t}`))return `<td class="outside-domain" aria-label="${t>n?'Not applicable: order exceeds the qubit number':'No catalogue record'}"><span aria-hidden="true">—</span></td>`;
      const entry=getEntry(n,t),s=status(entry),isSelected=n===selected.n&&t===selected.t;
      if(!matches(entry))return '<td class="filtered-cell" aria-label="Hidden by evidence filter">·</td>';
      matching++;
      return `<td><button class="state-cell ${s}${isSelected?' selected':''}" type="button" data-n="${n}" data-t="${t}" aria-pressed="${isSelected}" aria-label="${n} qubits, order ${t}: ${labels[s]}"><span class="cell-mark" aria-hidden="true">${symbols[s]}</span><span class="cell-label">${shortLabels[s]}</span></button></td>`;
    }).join('')}</tr>`).join('');
    $('state-table').innerHTML=`<table class="state-grid"><caption class="sr-only">Existence of pure anticoherent states by qubit number n and order t</caption>${head}<tbody>${body}</tbody></table>${matching ? '' : '<p class="no-matches">No cells match this evidence filter.</p>'}`;
    $('result-count').textContent=`${matching} ${matching===1?'cell':'cells'} · ${data.validation?.exactWitnessCount || states.size} exact representatives`;
    typeset($('state-table'));
  }
  function decimal(c) {
    const [re,im]=Array.isArray(c)?c:[Number(c),0];
    const f=x=>Math.abs(x)<1e-14?'0':Number(x).toPrecision(9).replace(/(?:\.0+|(\.\d*?)0+)$/,'$1');
    if (!im || Math.abs(im)<1e-14)return f(re);
    if (!re || Math.abs(re)<1e-14)return `${f(im)}i`;
    return `${f(re)} ${im<0?'−':'+'} ${f(Math.abs(im))}i`;
  }
  function exactText(expression) { return String(expression).replace(/\bI\b/g,'i').replace(/sqrt\(/g,'√(').replace(/\*\*/g,'^').replace(/\*/g,'·'); }
  function coefficients(state) {
    const support=state.support || (state.exactCoefficients||state.coefficients||[]).map((_,i)=>i).filter(i=>String(state.exactCoefficients?.[i])!=='0' && JSON.stringify(state.coefficients?.[i])!=='[0,0]');
    return `<details class="coefficient-disclosure"><summary>Dicke coefficients <span>${support.length} nonzero</span></summary><div class="coefficient-scroll"><table class="coefficient-table"><thead><tr><th scope="col">r</th><th scope="col">Exact amplitude</th><th scope="col">Decimal</th></tr></thead><tbody>${support.map(r=>`<tr><th scope="row">${r}</th><td class="exact-amplitude">${escape(state.exactCoefficients?.[r]!==undefined?exactText(state.exactCoefficients[r]):'—')}</td><td class="decimal-amplitude">${escape(state.coefficients?.[r]!==undefined?decimal(state.coefficients[r]):'—')}</td></tr>`).join('')}</tbody></table></div><p class="micro-copy">Amplitudes are ordered by excitation number. i denotes the imaginary unit; omitted amplitudes are zero.</p></details>`;
  }
  function renderDetail() {
    window.MajoranaViewer?.close();
    const entry=getEntry(selected.n,selected.t),s=status(entry);
    let candidates=(entry.stateIds||[]).map(id=>states.get(id)).filter(Boolean);
    let witness=candidates.find(x=>x.id===selected.stateId)||candidates[0];
    const majorana=witness && window.MAJORANA_DATA?.records?.[witness.id];
    selected.stateId=witness?.id||null;
    let content=`<div class="detail-heading"><p class="eyebrow">Selected pair</p>${badge(entry)}</div><h3>${selected.n} qubits <span class="detail-separator">/</span> order ${selected.t}</h3>`;
    if(s==='exists' && witness) {
      content+=`<div class="witness-metric">${math(`\\rho_{${selected.t}}=\\frac{I_{${selected.t+1}}}{${selected.t+1}}`,true)}</div>`;
      content+=`<label class="witness-label" for="witness-select">Representative <select id="witness-select">${candidates.map(st=>`<option value="${escape(st.id)}"${st.id===witness.id?' selected':''}>${escape((st.label||st.stateName||st.id)+symmetryLabel(st))}</option>`).join('')}</select></label>`;
      content+=symmetrySummary(witness);
      content+=`<p class="witness-order">${witness.achievedOrder}-AC representative<span>Witnesses orders 1–${witness.achievedOrder}</span></p>`;
      content+=`<div class="state-formula" role="region" tabindex="0" aria-label="Exact Dicke state formula">${math(witness.stateLatex||'',true)}</div>`;
      if(majorana?.points?.length)content+='<button id="view-majorana" class="view-majorana-button" type="button" aria-haspopup="dialog" aria-controls="majorana-viewer">View Majorana sphere ↗</button>';
      content+=`<div class="verification-line"><span class="verification-dot" aria-hidden="true"></span><span>Exact reduced-state verification${math(`\\quad\\operatorname{Tr}(\\rho_{${selected.t}}^2)=\\frac1{${selected.t+1}}`)}</span></div>`;
      if(witness.provenance)content+=`<p class="provenance">${escape(witness.provenance)}</p>`;
      content+=coefficients(witness);
      if(witness.notes)content+=`<details class="coefficient-details"><summary>Construction notes</summary><p class="detail-note">${escape(witness.notes)}</p></details>`;
    } else if(s==='impossible') {
      content+='<p class="result-statement">No pure symmetric state can be anticoherent to this order for these parameters.</p>';
      const proof=proofs.get(entry.proofId);
      if(proof)content+=`<div class="proof-summary"><h4>${escape(proof.title)}</h4><p>${escape(proof.summary||'')}</p>${proof.formulaLatex?math(proof.formulaLatex,true):''}<a class="text-link" href="#proof-${escape(proof.id)}">Read the argument <span aria-hidden="true">↘</span></a></div>`;
    } else if(s==='reported_absent') {
      content+='<p class="result-statement">Reported absent in the source.</p><div class="evidence-note">The catalogue records this published conclusion. A verified nonexistence proof is not supplied for this cell.</div>';
    } else if(s==='unknown') {
      content+='<p class="result-statement">Existence is unresolved in this catalogue.</p><div class="evidence-note">No verified witness or exclusion is recorded. This cell makes no claim based on an unsuccessful search.</div>';
    } else {
      content+='<p class="result-statement">The existence record is available; its representative data is not yet attached.</p>';
    }
    if(entry.notes && !['exists','impossible'].includes(s))content+=`<p class="detail-note">${escape(entry.notes)}</p>`;
    const refs=[...(entry.references||[]),...(witness?.references||[])];
    if(refs.length)content+=`<h4 class="detail-label">Source record</h4><ul class="detail-sources">${referenceLinks(refs)}</ul>`;
    content+='<div class="detail-actions"><button type="button" id="download-record">Download record ↓</button><button type="button" id="copy-link">Copy link ↗</button></div>';
    $('result-detail').innerHTML=content; typeset($('result-detail'));
    $('view-majorana')?.addEventListener('click',()=>window.MajoranaViewer?.open(majorana));
    $('witness-select')?.addEventListener('change',event=>{selected.stateId=event.target.value;renderDetail();updateHash();$('witness-select')?.focus({preventScroll:true});});
    $('download-record').addEventListener('click',()=>{
      const payload={entry,selectedState:witness||null,proof:proofs.get(entry.proofId)||null,majoranaRepresentation:majorana||null};
      const blob=new Blob([JSON.stringify(payload,null,2)+'\n'],{type:'application/json'});const url=URL.createObjectURL(blob);const a=document.createElement('a');a.href=url;a.download=`anticoherent-n${selected.n}-t${selected.t}.json`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);notify('Record downloaded');
    });
    $('copy-link').addEventListener('click',async()=>{
      const url=location.href.split('#')[0]+hashFor();
      try { if(!navigator.clipboard?.writeText)throw new Error('Clipboard unavailable');await navigator.clipboard.writeText(url);notify('Link copied'); }
      catch { const input=document.createElement('textarea');input.value=url;input.className='copy-fallback';document.body.append(input);input.select();const ok=document.execCommand('copy');input.remove();notify(ok?'Link copied':'Copy the address from your browser to share this result.'); }
    });
  }
  function renderProofs() {
    $('proof-list').innerHTML=(data.proofs||[]).map(p=>`<details class="proof-item" id="proof-${escape(p.id)}"><summary><span>${escape(p.title)}</span><span class="proof-expand" aria-hidden="true">+</span></summary><div class="proof-body"><p>${escape(p.summary||'')}</p>${p.formulaLatex?math(p.formulaLatex,true):''}${cleanBody(p.bodyHtml||'')}${p.references?.length?`<ul class="detail-sources">${referenceLinks(p.references)}</ul>`:''}</div></details>`).join('')||'<p class="muted">Proof records will appear here when included in the source data.</p>';
  }
  function renderConstructions() {
    const best=new Map();
    [...states.values()].forEach(st=>{if((st.n>=11 || st.achievedOrder>=4) && (!best.has(st.n)||st.achievedOrder>best.get(st.n).achievedOrder))best.set(st.n,st);});
    $('construction-list').innerHTML=[...best.values()].sort((a,b)=>a.n-b.n).map(st=>`<article class="construction-card"><p class="eyebrow">N = ${st.n} / order ${st.achievedOrder}</p><h3>${escape(st.label||st.stateName)}</h3>${symmetrySummary(st)}<p>${escape(st.provenance||'Exact representative with verified symmetric reductions.')}</p><button type="button" class="text-link" data-state-link="${escape(st.id)}">Inspect the state <span aria-hidden="true">↗</span></button></article>`).join('')||'<p class="muted">Source-linked constructions will appear here.</p>';
    $('general-constructions').innerHTML=(data.constructions||[]).map(c=>`<article class="general-construction"><h3>${escape(c.title)}</h3><p>${escape(c.text||'')}</p>${c.formulaLatex?`<div class="construction-formula">${math(c.formulaLatex,true)}</div>`:''}${c.references?.length?`<ul class="detail-sources">${referenceLinks(c.references)}</ul>`:''}</article>`).join('');
  }
  function renderSources() {
    const caveats=data.sourceCaveats||[];
    $('source-caveats').innerHTML=caveats.length?`<details class="caveat-box"><summary>Reading notes & source discrepancies <span>${caveats.length}</span></summary><div>${caveats.map(c=>`<p>${escape(typeof c==='string'?c:c.issue||c.notes||'')}${typeof c==='object'&&c.file?` <a href="${escape(sourceUrl(c))}">${escape(c.file)}${c.page?' · p. '+c.page:''} ↗</a>`:''}</p>`).join('')}</div></details>`:'';
    $('source-list').innerHTML=(data.sources||[]).map((src,i)=>`<article class="source-item"><span class="source-index">${String(i+1).padStart(2,'0')}</span><div><h3><a href="${escape(safeUrl(sourceUrl(src)))}">${escape(src.title||src.filename||src.file)} <span aria-hidden="true">↗</span></a></h3><p>${escape(Array.isArray(src.authors)?src.authors.join(', '):src.authors||'')}${src.year?' · '+escape(src.year):''}</p>${src.notes||src.description?`<p class="source-description">${escape(src.notes||src.description)}</p>`:''}<p class="source-filename">${escape(src.filename||src.file||'')}</p></div></article>`).join('');
    const definitionRefs=data.definitionReferences||[];
    $('definition-sources').innerHTML=definitionRefs.length?`Definition sources: ${definitionRefs.map(r=>`<a href="${escape(sourceUrl(r))}">${escape(r.file||r.title)}${r.page?' · p. '+r.page:''}</a>`).join('; ')}.`:'';
  }
  function handleAnchors() {
    if(location.hash.startsWith('#proof-')){ const proof=document.getElementById(location.hash.slice(1));if(proof?.tagName==='DETAILS')proof.open=true; }
  }
  $('state-table').addEventListener('click',event=>{const b=event.target.closest('button[data-n]');if(b)select(b.dataset.n,b.dataset.t,null,{focus:true});});
  $('state-table').addEventListener('keydown',event=>{
    const b=event.target.closest('button[data-n]');if(!b)return;
    const keys=['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','Home','End'];if(!keys.includes(event.key))return;
    event.preventDefault();const ns=numbers(),ts=orders();let ri=ns.indexOf(Number(b.dataset.n)),ci=ts.indexOf(Number(b.dataset.t));
    if(event.key==='Home')ci=0;else if(event.key==='End')ci=ts.length-1;else {ri+=event.key==='ArrowDown'?1:event.key==='ArrowUp'?-1:0;ci+=event.key==='ArrowRight'?1:event.key==='ArrowLeft'?-1:0;}
    if(ri<0||ri>=ns.length||ci<0||ci>=ts.length)return;
    const target=$('state-table').querySelector(`[data-n="${ns[ri]}"][data-t="${ts[ci]}"]`);target?.focus();
  });
  $('status-filter').addEventListener('change',event=>{filter=event.target.value;ensureSelection();renderMatrix();renderDetail();updateHash();});
  $('extended-toggle').addEventListener('change',event=>{extended=event.target.checked;ensureSelection();renderMatrix();renderDetail();updateHash();});
  $('construction-list').addEventListener('click',event=>{const b=event.target.closest('[data-state-link]');const st=states.get(b?.dataset.stateLink);if(st){select(st.n,st.achievedOrder,st.id);$('catalogue').scrollIntoView({behavior:'smooth',block:'start'});}});
  document.addEventListener('click',event=>{const link=event.target.closest('a[href^="#proof-"]');if(link){const p=document.getElementById(link.hash.slice(1));if(p)p.open=true;}});
  window.addEventListener('hashchange',()=>{if(location.hash.startsWith('#n=')){readHash();select(selected.n,selected.t,selected.stateId);}else handleAnchors();});
  if(typeof data.scope==='string')$('coverage-note').textContent=data.scope+' An unresolved cell makes no claim about existence.';
  readHash();syncControls();ensureSelection();renderMatrix();renderDetail();renderProofs();renderConstructions();renderSources();typeset();handleAnchors();
  // Keep direct section anchors intact on first load.
  if(!location.hash || location.hash.startsWith('#n='))updateHash();
})();
