/* Offline viewer. Scientific data and proof records live in data/catalogue.js. */
(() => {
  'use strict';
  const $ = id => document.getElementById(id);
  const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const math = (latex, display=false) => `<span class="math" data-tex="${escape(latex)}" data-display="${display}">${escape(latex)}</span>`;
  function scientific(value) {
    if(value==null)return 'Not supplied';
    const number=Number(value);if(!Number.isFinite(number))return escape(value);
    if(number===0)return math('0');
    const [mantissa,exponent]=number.toExponential(3).split('e');
    const coefficient=mantissa.replace(/(?:\.0+|(\.\d*?)0+)$/,'$1');
    return math(Number(exponent)===0?coefficient:`${coefficient}\\times10^{${Number(exponent)}}`);
  }
  function typeset(root=document) {
    root.querySelectorAll('[data-tex]').forEach(el => {
      if (window.katex) el.innerHTML=window.katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==='true',throwOnError:false,trust:false,strict:'ignore'});
    });
  }
  typeset();
  const data=window.MAJORANA_CATALOGUE;
  if (!data || !Array.isArray(data.entries)) {
    $('state-table').innerHTML='<p class="empty-state">The catalogue records could not be loaded. Keep the data folder alongside this page.</p>';
    $('result-count').textContent='Source records unavailable';
    return;
  }
  const numericalIds=new Set((data.numericalStates||[]).map(s=>s.id));
  const isNumerical=state=>Boolean(state&&(numericalIds.has(state.id)||state.numerical));
  const states=new Map([...(data.states||[]),...(data.numericalStates||[])].map(s=>[s.id,s]));
  const partitionKey=partition=>Array.isArray(partition)?[...partition].sort((a,b)=>b-a).join('-'):String(partition||'');
  const keyFor=entry=>`${entry.n}:${entry.partitionKey||partitionKey(entry.partition)}`;
  const entries=new Map(data.entries.map(e=>[keyFor(e),e]));
  const proofs=new Map((data.proofs||[]).map(p=>[p.id,p]));
  const sourceMap=new Map((data.sources||[]).map(s=>[s.filename,s]));
  const numbers=data.numbers||[2,3,4,5,6,7,8,9,10,11,12];
  $('n-filter').innerHTML='<option value="all">All</option>'+numbers.map(n=>`<option value="${Number(n)}">${Number(n)}</option>`).join('');
  const labels={proven:'Proven maximum',bounded:'Bounds only',no_ac:'No positive AC order'};
  let selected={n:12,partitionKey:Array(12).fill(1).join('-'),stateId:null};
  let filter='all',nFilter='12';
  let toastTimer;
  let disposeSphere=()=>{};
  let disposeViewerSphere=()=>{},viewerState=null,viewerCandidates=[],viewerReturnFocus=null;
  const symmetryCache=new Map();
  const groupCache=new Map();
  const symmetryColor=order=>({2:'#34786d',3:'#82609b',4:'#a45635',5:'#b38b23',continuous:'#356f94'}[order]||'#526b9b');

  function symmetryFor(state) {
    if(!symmetryCache.has(state.id)){
      let result;
      try{result=window.MajoranaSymmetry?.detect(state.majoranaPoints,{tolerance:1e-6});}
      catch{result=null;}
      symmetryCache.set(state.id,result||{status:'invalid',axes:[],notes:['Symmetry detection unavailable.']});
    }
    return symmetryCache.get(state.id);
  }
  function groupsFor(state) {
    if(!state?.majoranaPoints?.length)return null;
    if(!groupCache.has(state.id)){
      const unknown={symbol:null,name:'Not determined',kind:'unknown',order:null,numerical:true};
      const result=symmetryFor(state);
      const rotational=window.MajoranaGroups?.classify(result)||unknown;
      const point=window.MajoranaPointGroups?.classify(state.majoranaPoints,result,rotational)||unknown;
      groupCache.set(state.id,{point,rotational});
    }
    return groupCache.get(state.id);
  }
  const groupMath=group=>group?.symbol?math(group.symbol):'<span class="group-unknown">Not determined</span>';
  const groupPlain=group=>group?.symbol?group.symbol.replace(/\\infty/g,'∞').replace(/\\(?:mathrm|operatorname)/g,'').replace(/[{}_]/g,''):'unclassified';
  function groupCard(state,compact=false) {
    const groups=groupsFor(state);if(!groups)return '';
    const size=(group,unit)=>group.continuous?'Continuous':Number.isInteger(group.order)?`${group.order} ${unit}`:'Classification unavailable';
    return `<div class="constellation-groups${compact?' compact-groups':''}"><div class="group-item"><span class="group-kind">Point group</span><strong>${groupMath(groups.point)}</strong><span class="group-name">${escape(groups.point.name)}</span><span class="group-size">${escape(size(groups.point,'symmetries'))}</span></div><div class="group-item"><span class="group-kind">Rotational subgroup</span><strong>${groupMath(groups.rotational)}</strong><span class="group-name">${escape(groups.rotational.name)}</span><span class="group-size">${escape(size(groups.rotational,'rotations'))}</span></div><p class="group-evidence">Numerical geometry · This representative only. The point group includes any reflection or inversion symmetries.</p></div>`;
  }
  function stateOption(state,current) {
    return `<option value="${escape(state.id)}"${state.id===current?' selected':''}>${escape(state.label||state.id)} · ${isNumerical(state)?'reported order':'AC order'} ${state.achievedOrder} · ${escape(groupPlain(groupsFor(state)?.point))}</option>`;
  }
  function entryGroups(entry) {
    const grouped=new Map();
    for(const id of [...(entry.stateIds||[]),...(entry.numericalStateIds||[])]){
      const state=states.get(id);if(!state||partitionKey(state.partition)!==(entry.partitionKey||partitionKey(entry.partition)))continue;
      const group=groupsFor(state)?.point;if(group)grouped.set(group.symbol||'unknown',group);
    }
    if(!grouped.size)return '<span class="no-report" aria-label="No attached representative">—</span>';
    return `<div class="table-groups" title="Point groups of the attached representatives; the partition does not fix the group.">${[...grouped.values()].sort((a,b)=>(a.symbol||'').localeCompare(b.symbol||'')).map(group=>groupMath(group)).join('<span class="group-separator" aria-hidden="true">,</span>')}</div>`;
  }
  function prepareSymmetry(state) {
    const result=symmetryFor(state),counts=new Map();
    for(const axis of result.axes) {const key=axis.order==null?'continuous':String(axis.order);counts.set(key,(counts.get(key)||0)+1);}
    const keys=[...counts.keys()].sort((a,b)=>a==='continuous'?-1:b==='continuous'?1:Number(b)-Number(a));
    if(result.continuousFamilies?.length&&!keys.includes('2'))keys.push('2');
    $('viewer-symmetry-order').innerHTML='<option value="all">All orders</option>'+keys.map(key=>`<option value="${key}">${key==='continuous'?'Continuous axis':key+'-fold axes'}</option>`).join('');
    // Dense polyhedral groups are easier to read one rotation order at a time.
    $('viewer-symmetry-order').value=result.axes.length>8?keys[0]:'all';
    const available=Boolean(result.axes.length||result.continuousFamilies?.length);
    $('viewer-symmetry').checked=available;$('viewer-symmetry').disabled=!available;$('viewer-symmetry-order').disabled=!available;
    const family=result.continuousFamilies?.length?`<span class="symmetry-family">Every axis in the highlighted perpendicular plane is ${math('C_2')}.</span>`:'';
    $('viewer-symmetry-summary').innerHTML=available?keys.filter(key=>counts.has(key)).map(key=>`<span class="symmetry-key"><i style="background:${symmetryColor(key)}" aria-hidden="true"></i>${math('C_{'+(key==='continuous'?'\\infty':key)+'}')}<span>${counts.get(key)} ${counts.get(key)===1?'axis':'axes'}</span></span>`).join('')+family:result.status==='invalid'?'Symmetry axes could not be determined from these coordinates.':'No nontrivial rotational symmetry detected for this representative.';
    $('viewer-symmetry-method').innerHTML=`${math('C_m')} labels an axis for a rotation by ${math('2\\pi/m')} that preserves the whole weighted constellation. Each axis is labelled by its highest rotational order; opposite ends denote one axis. The rotational subgroup contains proper rotations only; the full point group also includes any reflections, inversion and rotation–reflection operations. Every operation must preserve point multiplicities. Groups are checked geometrically on the displayed coordinates, with distance tolerance ${math('10^{-6}')} on the unit sphere. This is a numerical symmetry check, separate from the AC certificates. ${result.status==='continuous'?'The continuous axis admits any rotation angle.':''}`;
    typeset($('viewer-symmetry-summary'));typeset($('viewer-symmetry-method'));
    return result;
  }

  function updateViewerAvailability(state,candidates=[]) {
    viewerState=state;viewerCandidates=candidates;
    const available=Boolean(state?.majoranaPoints?.length);
    $('view-majorana').disabled=!available;
    $('view-majorana').title=available?`${state.n} qubits · ${state.label||'selected state'}`:'No explicit state is attached to this partition.';
    $('majorana-availability').hidden=available;
    $('majorana-availability').textContent=available?'':'No explicit state attached for visualization.';
    if($('majorana-viewer').open){if(available)renderMajoranaViewer();else $('majorana-viewer').close();}
  }
  function renderMajoranaViewer() {
    const state=viewerState;if(!state?.majoranaPoints?.length)return;
    disposeViewerSphere();
    $('viewer-title').textContent=`${state.n} qubits · Majorana sphere`;
    $('viewer-summary').innerHTML=`${math(`\\lambda=${partitionLatex(state.partition)}`)}<span>${state.n} points · ${state.majoranaPoints.length} distinct positions</span>`;
    $('viewer-witness').innerHTML=viewerCandidates.map(s=>stateOption(s,state.id)).join('');
    $('viewer-groups').innerHTML=groupCard(state,true);typeset($('viewer-groups'));
    $('viewer-evidence').textContent=isNumerical(state)?`Reported AC order ${state.achievedOrder}, checked numerically. Coordinates visualize the refined saved state; its full multiplicity partition is certified.`:`Exact witness of AC order ${state.achievedOrder}. Coordinates are a numerical rendering of this representative.`;
    $('viewer-sphere').setAttribute('aria-label',`Majorana sphere: ${state.n} qubits, ${state.majoranaPoints.length} distinct positions. ${state.label||''}`);
    $('viewer-axes').checked=false;$('viewer-hull').checked=true;
    const symmetry=prepareSymmetry(state);
    typeset($('viewer-summary'));
    disposeViewerSphere=attachSphere($('viewer-sphere'),state.majoranaPoints,{
      download:$('viewer-download-svg'),exportLabel:state.label||state.id,reset:$('viewer-reset'),align:$('viewer-align'),hull:$('viewer-hull'),
      axes:$('viewer-axes'),
      symmetry,symmetryToggle:$('viewer-symmetry'),symmetryOrder:$('viewer-symmetry-order'),
      turns:[[$('viewer-left'),-.2,0],[$('viewer-right'),.2,0],[$('viewer-up'),0,-.2],[$('viewer-down'),0,.2]]
    });
  }
  function openMajoranaViewer() {
    if(!viewerState?.majoranaPoints?.length)return;
    viewerReturnFocus=document.activeElement;
    $('majorana-viewer').showModal();document.body.classList.add('majorana-viewer-open');
    renderMajoranaViewer();$('viewer-sphere').focus({preventScroll:true});
  }

  function safeUrl(value) {
    const url=String(value||'').trim();
    if (/^https?:\/\//i.test(url)) return url;
    return url && !/^[a-z][a-z\d+.-]*:/i.test(url) && !url.startsWith('//') ? url : '#sources';
  }
  function sourceUrl(ref) {
    if (ref.file || ref.filename) return `references/${encodeURIComponent(ref.file||ref.filename)}${ref.page?'#page='+Number(ref.page):''}`;
    if (ref.pdf) return safeUrl(ref.pdf.includes('/')?ref.pdf:'references/'+encodeURIComponent(ref.pdf));
    return safeUrl(ref.url||'#sources');
  }
  function referenceLinks(refs=[]) {
    const unique=new Map(refs.map(ref=>[`${ref.file||ref.filename||ref.url}:${ref.page||''}:${ref.location||''}`,ref]));
    return [...unique.values()].map(ref=>{
      const source=sourceMap.get(ref.file||ref.filename);
      return `<li><a class="source-link" href="${escape(sourceUrl(ref))}">${escape(ref.title||source?.title||ref.file||ref.filename||'Source')}${ref.page?' · p. '+escape(ref.page):''} ↗</a>${ref.location?`<span class="source-location">${escape(ref.location)}</span>`:''}</li>`;
    }).join('');
  }
  function cleanBody(html='') {
    const template=document.createElement('template');template.innerHTML=html;
    const allowed=new Set(['P','UL','OL','LI','EM','STRONG','CODE','PRE','SPAN','DIV','A','SUB','SUP','BR','DETAILS','SUMMARY','BLOCKQUOTE','H3','H4','H5','TABLE','THEAD','TBODY','TR','TH','TD','CAPTION']);
    [...template.content.querySelectorAll('*')].forEach(el=>{
      if (!allowed.has(el.tagName)) {el.replaceWith(document.createTextNode(el.textContent));return;}
      [...el.attributes].forEach(attr=>{
        if (!['class','data-tex','data-display','href','colspan','rowspan','scope'].includes(attr.name)) el.removeAttribute(attr.name);
        if (attr.name==='href') el.setAttribute('href',safeUrl(attr.value));
      });
    });
    return template.innerHTML;
  }
  function notify(message) { $('toast').textContent=message;$('toast').classList.add('visible');clearTimeout(toastTimer);toastTimer=setTimeout(()=>$('toast').classList.remove('visible'),2600); }
  function partitionLatex(partition) {
    const counts=new Map();[...(partition||[])].sort((a,b)=>b-a).forEach(m=>counts.set(m,(counts.get(m)||0)+1));
    return `(${[...counts].map(([m,count])=>count>1?`${m}^{${count}}`:String(m)).join(',')})`;
  }
  const entryFor=()=>entries.get(`${selected.n}:${selected.partitionKey}`);
  const matches=entry=>(filter==='all'||entry.status===filter)&&(nFilter==='all'||entry.n===Number(nFilter));
  const sortedEntries=[...data.entries].sort((a,b)=>a.n-b.n||b.lower-a.lower||b.upper-a.upper||(b.partition?.length||0)-(a.partition?.length||0)||keyFor(a).localeCompare(keyFor(b)));
  const visibleEntries=()=>sortedEntries.filter(matches);
  function hashFor() { const params=new URLSearchParams({n:selected.n,partition:selected.partitionKey});if(selected.stateId)params.set('state',selected.stateId);return '#'+params.toString(); }
  function updateHash() { if(!entryFor())return;const hash=hashFor();if(location.hash!==hash)history.replaceState(null,'',hash); }
  function readHash() {
    const params=new URLSearchParams(location.hash.slice(1)),n=Number(params.get('n'));
    if(!numbers.includes(n))return;
    let key=params.get('partition'),stateId=params.get('state');
    if(key)key=partitionKey(key.split('-').map(Number));
    if(!entries.has(`${n}:${key}`)&&params.has('t')){
      const legacy=(data.legacyEntries||[]).find(e=>e.n===n&&e.t===Number(params.get('t')));
      const legacyStates=(legacy?.stateIds||[]).map(id=>states.get(id)).filter(s=>s&&Array.isArray(s.partition)&&entries.has(`${n}:${partitionKey(s.partition)}`)).sort((a,b)=>b.achievedOrder-a.achievedOrder);
      const state=legacyStates.find(s=>s.id===stateId)||legacyStates[0];
      if(state){key=partitionKey(state.partition);stateId=state.id;}
    }
    if(!entries.has(`${n}:${key}`)){key=Array(n).fill(1).join('-');stateId=null;}
    if(entries.has(`${n}:${key}`)){selected={n,partitionKey:key,stateId};nFilter=String(n);}
  }
  function syncFilters(){$('n-filter').value=nFilter;$('status-filter').value=filter;}
  function ensureSelection() {
    if(!entryFor()||!matches(entryFor())){const first=visibleEntries()[0];if(first)selected={n:first.n,partitionKey:first.partitionKey||partitionKey(first.partition),stateId:null};}
  }
  function select(n,key,stateId=null,focus=false) {
    n=Number(n);if(!entries.has(`${n}:${key}`))return;
    selected={n,partitionKey:key,stateId};
    if(nFilter!=='all'&&nFilter!==String(n))nFilter=String(n);
    if(filter!=='all'&&entryFor().status!==filter)filter='all';
    syncFilters();renderTable();renderDetail();updateHash();
    if(focus)$('state-table').querySelector(`button[data-n="${n}"][data-partition="${key}"]`)?.focus({preventScroll:true});
  }
  function resultValue(entry) {
    if(entry.status==='no_ac')return '0';
    if(entry.status==='proven')return entry.value??entry.lower;
    return `${entry.lower}–${entry.upper}`;
  }
  function entryLabel(entry) {
    return `${entry.n} qubits, partition ${partitionLatex(entry.partition).replace(/\^\{(\d+)\}/g,' repeated $1')}: highest anticoherence order ${resultValue(entry)}, ${labels[entry.status]}${entry.reportedOrder!=null?', numerical source reports order '+entry.reportedOrder:''}`;
  }
  function renderTable() {
    const scrollTop=$('state-table').scrollTop,visible=visibleEntries();
    const rows=visible.map(entry=>{
      const key=entry.partitionKey||partitionKey(entry.partition),active=entry.n===selected.n&&key===selected.partitionKey;
      return `<tr class="partition-row ${escape(entry.status)}${active?' selected':''}" data-n="${entry.n}" data-partition="${escape(key)}"><td class="qubit-number">${entry.n}</td><th scope="row"><button type="button" class="partition-button" data-n="${entry.n}" data-partition="${escape(key)}" tabindex="${active?'0':'-1'}" aria-pressed="${active}" aria-label="${escape(entryLabel(entry))}">${math(entry.partitionLatex||partitionLatex(entry.partition))}</button></th><td class="order-value">${escape(resultValue(entry))}</td><td class="reported-value">${entry.reportedOrder!=null?`<span title="Numerical source report; not a rigorous bound">${escape(entry.reportedOrder)}</span>`:'<span class="no-report" aria-label="No numerical source report">—</span>'}</td><td class="group-cell">${entryGroups(entry)}</td><td><span class="table-evidence ${escape(entry.status)}"><i aria-hidden="true"></i>${entry.status==='proven'?'Proven':entry.status==='bounded'?'Bounds':'No positive AC'}</span></td></tr>`;
    }).join('');
    $('state-table').innerHTML=`<table class="partition-table"><caption class="sr-only">Highest anticoherence order by full multiplicity partition. Complete for 2–12 qubits; source-observed partitions at larger n. Numerical reports are separate from rigorous bounds. Point groups describe attached representatives, not entire partitions.</caption><thead><tr><th scope="col">${math('n')}</th><th scope="col">Partition ${math('\\lambda')}</th><th scope="col">${math('t_{\\max}')}</th><th scope="col" class="reported-heading">Reported<span>numerical</span></th><th scope="col" class="group-heading">Point groups<span>representatives</span></th><th scope="col">Evidence</th></tr></thead><tbody>${rows}</tbody></table>${visible.length?'':'<p class="no-matches">No partitions match these filters.</p>'}`;
    $('state-table').scrollTop=scrollTop;
    $('result-count').textContent=`${visible.length} of ${data.validation?.totalPartitionCount??data.validation?.partitionCount??data.entries.length} partitions · ${states.size} representatives`;
    typeset($('state-table'));
  }
  function decimal(pair) {
    const [re,im]=Array.isArray(pair)?pair:[Number(pair),0];
    const f=x=>Number(x)===0?'0':Number(x).toPrecision(9).replace(/(?:\.0+|(\.\d*?)0+)$/,'$1');
    if(!im)return f(re);
    if(!re)return f(im)+'i';
    return `${f(re)} ${im<0?'−':'+'} ${f(Math.abs(im))}i`;
  }
  function exactText(expression) { return String(expression).replace(/\bI\b/g,'i').replace(/sqrt\(/g,'√(').replace(/\*\*/g,'^').replace(/\*/g,'·'); }
  function coefficientTable(state) {
    const support=state.support||(state.coefficients||[]).map((c,r)=>[c,r]).filter(([c])=>Math.hypot(...c)>1e-13).map(([,r])=>r);
    return `<details class="coefficient-disclosure"><summary>Dicke coefficients <span>${support.length} ${isNumerical(state)?'amplitudes':'nonzero'}</span></summary><div class="coefficient-scroll"><table class="coefficient-table"><thead><tr><th scope="col">r</th><th scope="col">Exact amplitude</th><th scope="col">Decimal approximation</th></tr></thead><tbody>${support.map(r=>`<tr><th scope="row">${r}</th><td class="exact-amplitude">${state.coefficientLatex?.[r]!=null?math(state.coefficientLatex[r]):escape(state.exactCoefficients?.[r]!=null?exactText(state.exactCoefficients[r]):'—')}</td><td class="decimal-amplitude">${escape(state.coefficients?.[r]!=null?decimal(state.coefficients[r]):'—')}</td></tr>`).join('')}</tbody></table></div><p class="micro-copy">Normalized Dicke basis; r counts excitations. ${isNumerical(state)?'Displayed coefficients are numerical approximations.'+(state.verification?.exactMultiplicity&&(state.refinement||state.numerical?.refinement)?' The partition certificate applies to the full-precision saved vector.':''):'Omitted amplitudes are zero.'}${!isNumerical(state)&&!state.exactCoefficients&&state.verification?.certificatePath?' This state is specified by its linked certificate; the displayed decimals are approximations.':''}</p></details>`;
  }
  function pointTable(state) {
    if(!Array.isArray(state.majoranaPoints)||!state.majoranaPoints.length)return '';
    const f=x=>Number(x).toFixed(6).replace('-0.000000','0.000000');
    return `<details class="data-disclosure"><summary>Constellation coordinates <span>Numerical orientation</span></summary><table class="point-table"><thead><tr><th scope="col">Multiplicity</th><th scope="col">x</th><th scope="col">y</th><th scope="col">z</th></tr></thead><tbody>${state.majoranaPoints.map(p=>`<tr><td>${escape(p.multiplicity)}</td><td>${f(p.x)}</td><td>${f(p.y)}</td><td>${f(p.z)}</td></tr>`).join('')}</tbody></table></details>`;
  }
  function numericalChecks(state) {
    const report=state.numerical||{},check=report.ACCheck||{},cluster=report.clustering||{},criterion=report.criterion||{},refinement=state.refinement||report.refinement,certifiedPartition=Boolean(state.verification?.exactMultiplicity);
    const number=scientific;
    const stable=cluster.samePartitionAtTolerances||[];
    return `<details class="data-disclosure numerical-checks"${refinement?'':' open'}><summary>${refinement?'Original source checks':'Numerical verification'} <span>AC checked numerically</span></summary><dl><div><dt>Checked AC order</dt><dd>${escape(check.t??state.achievedOrder)}</dd></div><div><dt>Reduced-state Frobenius error</dt><dd>${number(check.frobeniusError)}</dd></div><div><dt>Squared Frobenius error</dt><dd>${number(check.purityExcess)}</dd></div><div><dt>Numerical acceptance threshold</dt><dd>${number(criterion.acceptanceThreshold)}</dd></div><div><dt>Point-clustering tolerance</dt><dd>${number(cluster.selectedTolerance)}</dd></div></dl><p class="micro-copy">AC errors compare the normalized reduced state with the maximally mixed symmetric state. ${report.claimNumericallySupported===false?'The reported order does not pass this numerical criterion.':'A small error supports numerical consistency, not an exact existence claim.'}</p><p class="micro-copy">${refinement?'Original source points were grouped':'Points are grouped'} by their Euclidean separation on the unit sphere. ${stable.length?`The same partition is recovered at ${stable.length} tested tolerances, from ${number(Math.min(...stable))} to ${number(Math.max(...stable))}.`:'The clustering-sensitivity record is linked below.'} ${certifiedPartition?'The full partition of the saved coefficients is independently certified by the linked root-count audit.':'This grouping does not certify exact repeated roots.'}</p></details>`;
  }
  function refinementDetails(state) {
    const refinement=state.refinement||state.numerical?.refinement;
    if(!refinement)return '';
    const number=scientific;
    return `<div class="refinement-note"><h4>Numerical accuracy refinement</h4><p class="refinement-residual"><span>${number(refinement.beforeResidual)}</span><span aria-label="after refinement">→</span><strong>${number(refinement.afterResidual)}</strong></p><p class="micro-copy">${escape(refinement.metric||'Residual before and after refinement at the reported AC order.')}${refinement.target!=null?` Target: ${number(refinement.target)}.`:''} Higher precision remains numerical evidence, not an exact AC certificate.</p>${refinement.partitionPreservation?`<p class="partition-certification">${escape(refinement.partitionPreservation)}</p>`:''}<div class="refinement-links">${refinement.coefficientPath?`<a class="text-link" href="${escape(safeUrl(refinement.coefficientPath))}">Refined coefficients at full precision ↗</a>`:''}${refinement.partitionCertificatePath?`<a class="text-link" href="${escape(safeUrl(refinement.partitionCertificatePath))}">Full-partition certificate ↗</a>`:''}</div></div>`;
  }
  function renderDetail() {
    disposeSphere();
    const entry=entryFor();
    if(!entry||!matches(entry)){
      updateViewerAvailability(null);
      $('result-detail').innerHTML='<p class="eyebrow">Selected partition</p><p class="empty-state">No partitions match these filters. Choose another qubit number or evidence category to explore the catalogue.</p>';
      return;
    }
    const key=entry.partitionKey||partitionKey(entry.partition);
    const candidates=[...(entry.stateIds||[]),...(entry.numericalStateIds||[])].map(id=>states.get(id)).filter(s=>s&&s.n===entry.n&&Array.isArray(s.partition)&&partitionKey(s.partition)===key).sort((a,b)=>Number(isNumerical(a))-Number(isNumerical(b))||b.achievedOrder-a.achievedOrder);
    const state=candidates.find(s=>s.id===selected.stateId)||candidates[0];
    selected.stateId=state?.id||null;
    const numerical=isNumerical(state);
    const certifiedPartition=Boolean(state?.verification?.exactMultiplicity);
    const refinement=state?.refinement||state?.numerical?.refinement;
    const entryProofs=[...new Set([...(entry.proofIds||[]),entry.proofId].filter(Boolean))].map(id=>proofs.get(id)).filter(Boolean);
    let html=`<div class="detail-heading"><p class="eyebrow">Selected partition</p><span class="badge ${escape(entry.status)}">${escape(labels[entry.status])}</span></div><h3>${entry.n} qubits</h3><div class="selected-partition">${math(`\\lambda=${entry.partitionLatex||partitionLatex(entry.partition)}`)}<span>${entry.partition.length} distinct point${entry.partition.length===1?'':'s'}</span></div>`;
    html+=`<div class="maximum-metric"><span class="maximum-value">${escape(resultValue(entry))}</span><div class="maximum-caption"><strong>${entry.status==='bounded'?'Bounds on the highest AC order':'Highest achievable AC order'}</strong><p>${entry.status==='no_ac'?'No positive AC order.':entry.status==='proven'?'Existence and upper bound agree.':'The exact maximum remains undetermined.'}</p></div></div>`;
    if(entry.coverage==='source-library')html+='<p class="coverage-detail">Source-library coverage: only observed multiplicity patterns are listed at this qubit number. This is not an enumeration of all partitions.</p>';
    if(entry.reportedOrder!=null)html+=`<div class="reported-order-note"><span>Numerical source report</span><strong>Order ${escape(entry.reportedOrder)}</strong><p>A reported candidate order, separate from the rigorous maximum or interval above.</p></div>`;
    if(entry.status==='no_ac')html+='<p class="detail-note">States with this partition exist, but none is 1-anticoherent. By convention, their highest achievable AC order is 0.</p>';
    if(state) {
      html+=`<label class="witness-label" for="witness-select">Representative for this partition<select id="witness-select">${candidates.map(s=>stateOption(s,state.id)).join('')}</select></label><div class="state-meta"><span>${numerical?'Reported numerical order':'Witness order'}: ${state.achievedOrder}</span><span>Largest multiplicity: ${state.maxMultiplicity}</span>${state.distinctPoints!=null?`<span>${state.distinctPoints} distinct points</span>`:''}</div>`;
      html+=groupCard(state);
      if(state.majoranaPoints?.length)html+='<button id="view-majorana-detail" class="view-majorana-button detail-viewer-button" type="button" aria-haspopup="dialog" aria-controls="majorana-viewer">View Majorana sphere ↗</button>';
      if(numerical)html+=`<p class="partition-note numerical-note">${certifiedPartition?'This representative has a certified full partition for its saved coefficients. Its AC order is checked numerically; this refinement does not prove exact anticoherence or optimality.':'This numerical representative illustrates the source report. It does not certify exact anticoherence, exact point multiplicities or attainment of a rigorous bound.'}</p>`;
      else if(entry.status==='proven')html+=`<p class="partition-note">${state.achievedOrder===(entry.value??entry.lower)?'This representative attains the maximum AC order for this partition.':`This representative has AC order ${state.achievedOrder}; the maximum for this partition is ${entry.value??entry.lower}.`}</p>`;
      else if(entry.status==='bounded')html+=`<p class="partition-note">${state.achievedOrder===entry.lower?'This representative attains the stated lower bound.':`This representative has AC order ${state.achievedOrder}; the established lower bound is ${entry.lower}.`} The upper bound is ${entry.upper}.</p>`;
      if(state.stateLatex)html+=`<div class="state-formula" role="region" tabindex="0" aria-label="Dicke state formula">${math(state.stateLatex,true)}</div>`;
      html+=`<div class="partition-line"><span>${certifiedPartition?'Certified full partition':numerical?'Numerically grouped partition':'Full partition'}</span>${math(`\\lambda=${partitionLatex(state.partition)}`)}</div>`;
      if(Array.isArray(state.majoranaPoints)&&state.majoranaPoints.length) {
        html+=`<div class="sphere-frame"><canvas id="majorana-sphere" class="sphere-canvas" tabindex="0" aria-label="Rotatable Majorana sphere for ${escape(state.label||'the selected state')}. Numbers label point multiplicities. Drag or use arrow keys to rotate; Home resets."></canvas><div class="sphere-controls"><span>Drag to rotate · Arrow keys when focused</span><label><input type="checkbox" id="inline-hull"> Convex hull</label><button id="inline-download-svg" type="button">Download SVG</button><button id="reset-sphere" type="button">Reset view</button></div><div class="sphere-legend"><span><i class="sphere-dot" aria-hidden="true"></i>Front hemisphere</span><span><i class="sphere-dot back" aria-hidden="true"></i>Back hemisphere</span><span><b class="mv-centre-key" aria-hidden="true">+</b>Sphere centre</span><span><b class="mv-barycentre-key" aria-hidden="true">◇</b>Barycentre (with multiplicities)</span><span>${state.majoranaPoints.length>30?'Unlabelled = 1':'Numbers beside points = multiplicities'}</span></div></div><p class="sphere-help">A representative orientation. Repeated points are drawn as one labelled cluster.${state.majoranaPoints.length>30?' Single points are smaller and unlabelled in this dense constellation.':''}</p>`;
      } else html+='<p class="detail-note">Constellation coordinates are not attached to this representative.</p>';
      const verification=state.verification||{};
      html+=`<div class="verification-line${verification.exactAC?'':' numerical'}"><span class="verification-dot" aria-hidden="true"></span><span>${verification.exactAC?'Anticoherence verified exactly':state.numerical?.claimNumericallySupported===false?'Reported AC order fails this numerical check':'AC checked numerically'}${verification.exactMultiplicity?numerical?' · Full partition certified for saved coefficients':' · Full partition certified':''}</span></div>`;
      if(verification.certificatePath)html+=`<p class="micro-copy"><a href="${escape(safeUrl(verification.certificatePath))}">${numerical?refinement?'Original numerical check record':'Numerical check record':'Verification certificate'} ↗</a>${verification.acCertificatePath?` · <a href="${escape(safeUrl(verification.acCertificatePath))}">AC verification ↗</a>`:''}</p>`;
      if(numerical)html+=numericalChecks(state)+refinementDetails(state);
      html+=coefficientTable(state);
      if(state.sourceFiles?.length)html+=`<details class="data-disclosure"><summary>Source data files</summary><ul class="raw-source-files">${state.sourceFiles.map(file=>`<li><a href="${escape(safeUrl(file.href))}">${escape(file.label||'Source data')} ↗</a></li>`).join('')}</ul><p class="micro-copy">Local counterparts matched to the saved primary-page links. Online byte identity could not be checked while the server was unreachable.</p></details>`;
      if(state.polynomialLatex)html+=`<details class="data-disclosure"><summary>Majorana polynomial</summary>${math(state.polynomialLatex,true)}<p class="micro-copy">An overall nonzero scalar does not change the roots.</p></details>`;
      if(state.algebraicDefinition){const a=state.algebraicDefinition;html+=`<details class="data-disclosure"><summary>Exact algebraic parameters</summary><p class="micro-copy">The parameter ${math('k')} is the unique real root in the interval below. These equations define the exact state; the displayed decimals approximate it.</p>${math(a.polynomialLatex,true)}${math(a.intervalLatex,true)}${math(a.recoveryLatex,true)}</details>`;}
      html+=pointTable(state);
      if(state.notes||verification.method)html+=`<details class="data-disclosure"><summary>Construction and verification</summary>${state.notes?`<p class="detail-note">${escape(state.notes)}</p>`:''}${verification.method?`<p class="detail-note">${escape(verification.method)}</p>`:''}</details>`;
    }
    if(entry.lower>0&&(!state||(numerical&&entry.attainment==='existence-theorem'))) {
      html+=`<div class="existence-note"><h4>Existence established without an attached vector</h4><p>${entry.attainment==='existence-theorem'?`A theorem guarantees a state with this exact partition and AC order at least ${entry.lower}.`:`The lower bound establishes AC order at least ${entry.lower} within this exact partition.`} An exact representative is not included in this record. The supporting existence and upper-bound arguments are linked below.</p></div>`;
    }
    if(entry.boundLatex)html+=`<div class="partition-bound">${math(entry.boundLatex,true)}</div>`;
    if(entryProofs.length)html+=`<div class="evidence-box"><h4>Evidence for this partition</h4>${entryProofs.map(proof=>`<div class="evidence-item"><a class="text-link" href="#proof-${escape(proof.id)}">${escape(proof.title)} ↘</a>${proof.summary?`<p>${escape(proof.summary)}</p>`:''}</div>`).join('')}</div>`;
    if(entry.upperReasons?.length)html+=`<details class="data-disclosure"><summary>Excluded higher orders</summary><ul class="upper-reasons">${entry.upperReasons.map(item=>`<li><strong>Order ${escape(item.order)}:</strong> ${escape(item.reason||'Excluded by the linked bound.')}${proofs.has(item.proofId)?` <a href="#proof-${escape(item.proofId)}">Proof ↘</a>`:''}</li>`).join('')}</ul></details>`;
    if(entry.notes)html+=`<p class="detail-note">${escape(entry.notes)}</p>`;
    const refs=[...(entry.references||[]),...(state?.references||[])];
    if(refs.length)html+=`<h4 class="detail-label">Source record</h4><ul class="detail-sources">${referenceLinks(refs)}</ul>`;
    html+='<div class="detail-actions"><button id="download-record" type="button">Download record ↓</button><button id="copy-link" type="button">Copy link ↗</button></div>';
    $('result-detail').innerHTML=html;typeset($('result-detail'));
    if(state&&$('majorana-sphere'))disposeSphere=attachSphere($('majorana-sphere'),state.majoranaPoints,{download:$('inline-download-svg'),exportLabel:state.label||state.id,reset:$('reset-sphere'),hull:$('inline-hull'),symmetry:symmetryFor(state)});
    updateViewerAvailability(state,candidates);
    $('view-majorana-detail')?.addEventListener('click',openMajoranaViewer);
    $('witness-select')?.addEventListener('change',event=>{selected.stateId=event.target.value;renderDetail();updateHash();$('witness-select')?.focus({preventScroll:true});});
    $('download-record').addEventListener('click',()=>{
      const blob=new Blob([JSON.stringify({entry,selectedState:state||null,symmetryGroups:state?groupsFor(state):null,proofs:entryProofs},null,2)+'\n'],{type:'application/json'});
      const url=URL.createObjectURL(blob);const a=document.createElement('a');a.href=url;a.download=`majorana-n${entry.n}-partition-${key}.json`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);notify('Record downloaded');
    });
    $('copy-link').addEventListener('click',async()=>{
      const url=location.href.split('#')[0]+hashFor();
      try{if(!navigator.clipboard?.writeText)throw new Error();await navigator.clipboard.writeText(url);notify('Link copied');}
      catch{const input=document.createElement('textarea');input.value=url;input.className='copy-fallback';document.body.append(input);input.select();const ok=document.execCommand('copy');input.remove();notify(ok?'Link copied':'Copy the address from your browser to share this result.');}
    });
  }

  // Axes are unoriented. Prefer the northern end, then y/x at the equator.
  // Equal orders retain the detector's stable ordering; a continuous axis wins.
  function principalAxis(symmetry){
    let chosen=null,best=1;
    for(const a of symmetry?.axes||[]){
      const rank=a.order==null?Infinity:Number(a.order),v=a.direction;
      if(!Array.isArray(v)||v.length!==3||!v.every(Number.isFinite)||Math.hypot(...v)<1e-12||rank<=best)continue;
      const norm=Math.hypot(...v),direction=v.map(x=>x/norm);
      const sign=[direction[2],direction[1],direction[0]].find(x=>Math.abs(x)>1e-12)<0?-1:1;
      chosen={order:a.order,direction:direction.map(x=>sign*x)};best=rank;
    }
    return chosen;
  }
  function axisToZ(direction){
    const norm=Math.hypot(...direction),[x,y,z]=direction.map(v=>v/norm);
    // Shortest proper rotation from direction to +z; handle the antipode.
    if(z<-1+1e-12)return [0,1,0,0];
    const q=[1+z,y,-x,0],length=Math.hypot(...q);
    return q.map(v=>v/length);
  }
  function attachSphere(canvas,points,controls={}) {
    const canvasContext=canvas.getContext('2d');if(!canvasContext)return ()=>{};
    const dense=points.length>30;
    const weight=points.reduce((sum,p)=>sum+(p.multiplicity??1),0);
    const barycentre=Object.fromEntries(['x','y','z'].map(axis=>[axis,points.reduce((sum,p)=>sum+(p.multiplicity??1)*p[axis],0)/weight]));
    const hull=window.MajoranaHull?.build(points);
    if(controls.hull){controls.hull.disabled=!hull?.faces.length;controls.hull.checked=Boolean(hull?.faces.length);controls.hull.title=hull?.faces.length?'Convex hull of the distinct Majorana positions':'At least three non-collinear positions are needed';}
    const multiply=(a,b)=>[a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0]];
    const axis=(x,y,z,angle)=>[Math.cos(angle/2),x*Math.sin(angle/2),y*Math.sin(angle/2),z*Math.sin(angle/2)];
    const principal=principalAxis(controls.symmetry);
    const initial=principal?axis(1,0,0,-Math.PI/2+.28):multiply(axis(1,0,0,-.28),axis(0,1,0,.42));
    const initialFrame=principal?axisToZ(principal.direction):[1,0,0,0];
    // Choose one end in the aligned state frame, independent of the camera.
    const axisDirections=new Map();
    for(const a of controls.symmetry?.axes||[]){
      const v=multiply(multiply(initialFrame,[0,...a.direction]),[initialFrame[0],-initialFrame[1],-initialFrame[2],-initialFrame[3]]);
      const sign=[v[3],v[2],v[1]].find(x=>Math.abs(x)>1e-12)<0?-1:1;
      axisDirections.set(a.direction,a.direction.map(x=>sign*x));
    }
    let orientation=[...initial],coordinateFrame=[...initialFrame],drag=null,frame=0,disposed=false;
    const listeners=[];
    const listen=(element,event,handler)=>{if(element){element.addEventListener(event,handler);listeners.push(()=>element.removeEventListener(event,handler));}};
    const rotate=(point,inFrame=true)=>{
      const q=inFrame?multiply(orientation,coordinateFrame):orientation;
      const v=multiply(multiply(q,[0,point.x,point.y,point.z]),[q[0],-q[1],-q[2],-q[3]]);
      return {x:v[1],y:v[2],z:v[3],multiplicity:point.multiplicity};
    };
    function turn(dx,dy){orientation=multiply(multiply(axis(0,1,0,dx),axis(1,0,0,dy)),orientation);const norm=Math.hypot(...orientation);orientation=orientation.map(v=>v/norm);schedule();}
    function schedule(){if(!frame&&!disposed)frame=requestAnimationFrame(()=>{frame=0;draw();});}
    function draw(context=canvasContext,exportBox=null){
      const box=exportBox||canvas.getBoundingClientRect();if(!box.width||!box.height)return;
      const dpr=exportBox?1:Math.min(window.devicePixelRatio||1,2);if(!exportBox){canvas.width=Math.round(box.width*dpr);canvas.height=Math.round(box.height*dpr);}context.setTransform(dpr,0,0,dpr,0,0);
      const width=box.width,height=box.height,cx=width/2,cy=height/2,radius=Math.min(width*.365,height*.405)*.95;
      context.clearRect(0,0,width,height);
      context.save();context.translate(cx,cy+radius*.95);context.scale(1,.12);const shadow=context.createRadialGradient(0,0,1,0,0,radius);shadow.addColorStop(0,'rgba(83,61,40,.12)');shadow.addColorStop(1,'rgba(83,61,40,0)');context.fillStyle=shadow;context.beginPath();context.arc(0,0,radius,0,Math.PI*2);context.fill();context.restore();
      const fill=context.createRadialGradient(cx-radius*.32,cy-radius*.38,radius*.06,cx,cy,radius);fill.addColorStop(0,'#fffefa');fill.addColorStop(.75,'#f3eee4');fill.addColorStop(1,'#e8dfd1');context.fillStyle=fill;context.beginPath();context.arc(cx,cy,radius,0,Math.PI*2);context.fill();
      const project=p=>({x:cx+radius*p.x,y:cy-radius*p.y,z:p.z});
      const circles=[];
      for(const latitude of [-Math.PI/3,-Math.PI/6,0,Math.PI/6,Math.PI/3])circles.push(Array.from({length:181},(_,i)=>{const a=i*Math.PI/90;return rotate({x:Math.cos(latitude)*Math.cos(a),y:Math.sin(latitude),z:Math.cos(latitude)*Math.sin(a)},false);}));
      for(const longitude of [0,Math.PI/4,Math.PI/2,3*Math.PI/4])circles.push(Array.from({length:181},(_,i)=>{const a=i*Math.PI/90;return rotate({x:Math.cos(a)*Math.cos(longitude),y:Math.sin(a),z:Math.cos(a)*Math.sin(longitude)},false);}));
      for(const front of [false,true]){
        context.strokeStyle=front?'rgba(133,112,88,.28)':'rgba(133,112,88,.16)';context.lineWidth=.8;context.setLineDash(front?[]:[2.5,4]);
        for(const circle of circles){context.beginPath();for(let i=1;i<circle.length;i++){const a=circle[i-1],b=circle[i];if((.5*(a.z+b.z)>=0)!==front)continue;const pa=project(a),pb=project(b);context.moveTo(pa.x,pa.y);context.lineTo(pb.x,pb.y);}context.stroke();}
      }
      context.setLineDash([]);context.strokeStyle='#c8b9a6';context.lineWidth=1;context.beginPath();context.arc(cx,cy,radius,0,Math.PI*2);context.stroke();
      if(controls.hull?.checked&&hull?.faces.length){
        const vertices=hull.vertices.map(point=>rotate(point));
        const faces=hull.faces.map(ids=>({ids,depth:ids.reduce((s,i)=>s+vertices[i].z,0)/ids.length})).sort((a,b)=>a.depth-b.depth);
        context.save();context.setLineDash([]);
        for(const face of faces){
          context.fillStyle=face.depth>=0?'rgba(58,126,119,.16)':'rgba(58,126,119,.07)';
          context.beginPath();face.ids.forEach((i,j)=>{const p=project(vertices[i]);if(j)context.lineTo(p.x,p.y);else context.moveTo(p.x,p.y);});context.closePath();context.fill();
        }
        context.lineWidth=1.1;
        for(const [a,b] of hull.edges){
          context.strokeStyle=(vertices[a].z+vertices[b].z)/2>=0?'rgba(45,111,105,.7)':'rgba(45,111,105,.28)';
          const p=project(vertices[a]),q=project(vertices[b]);context.beginPath();context.moveTo(p.x,p.y);context.lineTo(q.x,q.y);context.stroke();
        }
        context.restore();
      }
      const axisLabels=[];
      function axisLine(vector,label,color,unoriented=true){
        // Keep the saved end of an axis: changing its sign with the viewing
        // depth makes the label jump to the opposite end during rotation.
        const direction=unoriented?(axisDirections.get(vector)||vector):vector;
        const endpoint=rotate({x:direction[0],y:direction[1],z:direction[2]},unoriented);
        const p=project(endpoint),dx=p.x-cx,dy=p.y-cy,length=Math.hypot(dx,dy);
        // Extend in 3D, then project. An axis pointing towards the camera
        // must foreshorten with the constellation, not span a full diameter.
        const scale=unoriented?1.08:1.05,extent=length*scale;
        context.save();context.strokeStyle=color;context.lineWidth=unoriented?1.8:1;context.globalAlpha=.8;
        for(const sign of [-1,1]){context.setLineDash(sign*endpoint.z<0?[4,4]:[]);context.beginPath();context.moveTo(cx,cy);context.lineTo(cx+sign*dx*scale,cy+sign*dy*scale);context.stroke();}
        const endOn=unoriented?Math.max(0,1-length/10):0;
        if(endOn>0){context.globalAlpha=.8*endOn;context.setLineDash([]);context.beginPath();context.arc(cx,cy,4,0,2*Math.PI);context.stroke();}
        context.restore();
        // A continuous screen-space offset stays readable near an end-on
        // axis without dividing by its almost-zero projected length.
        const padding=unoriented?12:28,offset=padding/(extent+padding);
        axisLabels.push({label,color,x:Math.max(24,Math.min(width-24,cx+dx*scale*(1+offset))),y:Math.max(13,Math.min(height-13,cy+dy*scale*(1+offset)-padding*offset))});
      }
      if(controls.axes?.checked)for(const [name,vector,color] of [['x',[1,0,0],'#a15435'],['y',[0,1,0],'#41776e'],['z',[0,0,1],'#76648a']])axisLine(vector,name,color,false);
      if(controls.symmetryToggle?controls.symmetryToggle.checked:Boolean(controls.symmetry?.axes?.length)){
        const order=controls.symmetryOrder?.value||'all';
        for(const a of controls.symmetry.axes){const key=a.order==null?'continuous':String(a.order);if(order!=='all'&&order!==key)continue;const label=key==='continuous'?'C∞':'C'+key.replace(/\d/g,n=>'₀₁₂₃₄₅₆₇₈₉'[Number(n)]);axisLine(a.direction,label,symmetryColor(key));}
        if(controls.symmetry.continuousFamilies?.length&&(order==='all'||order==='2')){
          const normal=controls.symmetry.axes.find(a=>a.order==null)?.direction;
          if(normal){
            const helper=Math.abs(normal[0])<.8?[1,0,0]:[0,1,0];
            const u=[normal[1]*helper[2]-normal[2]*helper[1],normal[2]*helper[0]-normal[0]*helper[2],normal[0]*helper[1]-normal[1]*helper[0]],length=Math.hypot(...u);for(let i=0;i<3;i++)u[i]/=length;
            const v=[normal[1]*u[2]-normal[2]*u[1],normal[2]*u[0]-normal[0]*u[2],normal[0]*u[1]-normal[1]*u[0]];
            context.save();context.strokeStyle=symmetryColor('2');context.lineWidth=2.2;
            for(let i=0;i<180;i++){const point=angle=>rotate({x:u[0]*Math.cos(angle)+v[0]*Math.sin(angle),y:u[1]*Math.cos(angle)+v[1]*Math.sin(angle),z:u[2]*Math.cos(angle)+v[2]*Math.sin(angle)}),a=point(i*Math.PI/90),b=point((i+1)*Math.PI/90),pa=project(a),pb=project(b);context.setLineDash(a.z<0?[3,4]:[]);context.beginPath();context.moveTo(pa.x,pa.y);context.lineTo(pb.x,pb.y);context.stroke();}context.restore();
          }
        }
      }
      const ordered=points.map(point=>rotate(point)).sort((a,b)=>a.z-b.z);
      for(const point of ordered){
        const p=project(point),front=point.z>=0,labelled=Number(point.multiplicity)>1,nodeRadius=4.2+(point.z+1)*.6;
        context.save();context.globalAlpha=front?1:.85;context.lineWidth=front?1.3:1.4;context.setLineDash([]);
        context.fillStyle=front?'#a45a34':'#ddb392';context.strokeStyle=front?'#814321':'#8a4d2c';
        if(front){context.shadowColor='rgba(80,43,23,.19)';context.shadowBlur=3;context.shadowOffsetY=1;}
        context.beginPath();context.arc(p.x,p.y,nodeRadius,0,Math.PI*2);context.fill();context.stroke();context.shadowBlur=0;context.shadowOffsetY=0;context.setLineDash([]);
        if(labelled){const text=String(point.multiplicity),x=p.x+nodeRadius+5,y=p.y-nodeRadius-3;context.textAlign='left';context.textBaseline='middle';context.font='600 12px system-ui,-apple-system,sans-serif';context.lineJoin='round';context.lineWidth=2.5;context.strokeStyle='rgba(255,253,247,.85)';context.strokeText(text,x,y);context.fillStyle=front?'#703b20':'#744b33';context.fillText(text,x,y);}context.restore();
      }
      for(const label of axisLabels){context.save();context.font='600 13px system-ui,-apple-system,sans-serif';context.textAlign='center';context.textBaseline='middle';context.lineJoin='round';context.lineWidth=3;context.strokeStyle='rgba(255,253,247,.82)';context.strokeText(label.label,label.x,label.y);context.fillStyle=label.color;context.fillText(label.label,label.x,label.y);context.restore();}
      // Draw both markers on top of the translucent geometry. The cross
      // remains visible inside the diamond when their projections coincide.
      context.save();context.globalAlpha=1;context.shadowBlur=0;context.setLineDash([]);context.lineJoin='round';
      context.beginPath();context.moveTo(cx-4,cy);context.lineTo(cx+4,cy);context.moveTo(cx,cy-4);context.lineTo(cx,cy+4);
      context.strokeStyle='#fffdf7';context.lineWidth=4;context.stroke();context.strokeStyle='#35424a';context.lineWidth=1.7;context.stroke();
      const bary=project(rotate(barycentre));
      context.beginPath();context.moveTo(bary.x,bary.y-6);context.lineTo(bary.x+6,bary.y);context.lineTo(bary.x,bary.y+6);context.lineTo(bary.x-6,bary.y);context.closePath();
      context.strokeStyle='#fffdf7';context.lineWidth=4;context.stroke();context.strokeStyle='#a12b75';context.lineWidth=1.8;context.stroke();context.restore();
    }
    const down=event=>{if(event.button!==0&&event.pointerType==='mouse')return;canvas.focus({preventScroll:true});drag={id:event.pointerId,x:event.clientX,y:event.clientY};canvas.setPointerCapture(event.pointerId);canvas.classList.add('dragging');};
    const move=event=>{if(!drag||drag.id!==event.pointerId)return;const dx=event.clientX-drag.x,dy=event.clientY-drag.y;drag.x=event.clientX;drag.y=event.clientY;turn(dx*.008,dy*.008);};
    const up=event=>{if(drag&&drag.id===event.pointerId){drag=null;canvas.classList.remove('dragging');if(canvas.hasPointerCapture(event.pointerId))canvas.releasePointerCapture(event.pointerId);}};
    const reset=()=>{orientation=[...initial];coordinateFrame=[...initialFrame];schedule();};
    if(controls.align){
      controls.align.disabled=!principal;
      controls.align.title=principal?`Align ${principal.order==null?'the continuous':`a ${principal.order}-fold`} symmetry axis with z in the displayed frame`:'No rotational symmetry axis detected';
    }
    listen(controls.align,'click',()=>{
      if(!principal)return;
      coordinateFrame=axisToZ(principal.direction);
      // z is almost vertical, slightly tilted towards the viewer for depth.
      orientation=axis(1,0,0,-Math.PI/2+.28);
      if(controls.symmetryToggle)controls.symmetryToggle.checked=true;
      if(controls.symmetryOrder)controls.symmetryOrder.value=principal.order==null?'continuous':String(principal.order);
      schedule();
    });
    if(controls.download)controls.download.disabled=!window.MajoranaSVG;
    listen(controls.download,'click',()=>{
      const box=canvas.getBoundingClientRect();if(!box.width||!box.height)return;
      const svg=new window.MajoranaSVG.Context();draw(svg,box);
      const n=points.reduce((sum,p)=>sum+p.multiplicity,0);
      const label=controls.exportLabel||'Majorana constellation';
      const filename=`majorana-n${n}-${label.normalize('NFKD').replace(/[^a-zA-Z0-9]+/g,'-').replace(/^-|-$/g,'').slice(0,80)||'state'}.svg`;
      window.MajoranaSVG.download(svg.toSVG(box.width,box.height,`${n} qubits · ${label}`,'Majorana constellation in the current viewing orientation. Visible hull and axes follow the viewer settings. Numbers denote multiplicities greater than one. The dark cross marks the sphere centre; the magenta diamond marks the barycentre, averaging all Majorana vectors with their multiplicities.'),filename);
    });
    const key=event=>{const step=event.shiftKey ? 0.25 : 0.08;const directions={ArrowLeft:[-step,0],ArrowRight:[step,0],ArrowUp:[0,-step],ArrowDown:[0,step]};if(event.key==='Home'){event.preventDefault();reset();}else if(directions[event.key]){event.preventDefault();turn(...directions[event.key]);}};
    listen(canvas,'pointerdown',down);listen(canvas,'pointermove',move);listen(canvas,'pointerup',up);listen(canvas,'pointercancel',up);listen(canvas,'lostpointercapture',up);listen(canvas,'keydown',key);listen(controls.reset,'click',reset);
    listen(controls.hull,'change',schedule);listen(controls.axes,'change',schedule);
    listen(controls.symmetryToggle,'change',schedule);listen(controls.symmetryOrder,'change',schedule);
    for(const [button,dx,dy] of controls.turns||[])listen(button,'click',()=>turn(dx,dy));
    const observer=typeof ResizeObserver!=='undefined'?new ResizeObserver(schedule):null;observer?.observe(canvas);window.addEventListener('resize',schedule);schedule();
    return ()=>{disposed=true;observer?.disconnect();window.removeEventListener('resize',schedule);if(frame)cancelAnimationFrame(frame);listeners.forEach(remove=>remove());if(drag&&canvas.hasPointerCapture(drag.id))canvas.releasePointerCapture(drag.id);drag=null;canvas.classList.remove('dragging');};
  }

  function renderProofs(){
    $('proof-list').innerHTML=(data.proofs||[]).map(proof=>`<details class="proof-item" id="proof-${escape(proof.id)}"><summary>${escape(proof.title)}</summary><div class="proof-body">${proof.summary?`<p>${escape(proof.summary)}</p>`:''}${proof.formulaLatex?math(proof.formulaLatex,true):''}${cleanBody(proof.bodyHtml||'')}${proof.references?.length?`<ul class="detail-sources">${referenceLinks(proof.references)}</ul>`:''}</div></details>`).join('')||'<p class="section-intro">Proof records are not attached to this dataset.</p>';
  }
  function renderLibraryScope(){
    const library=data.library;if(!library){$('library-scope').hidden=true;return;}
    const upper=Math.max(...numbers);
    $('library-scope').hidden=false;
    $('library-scope').innerHTML=`<strong>Complete core: all 270 partitions for 2–12 qubits.</strong><p>For 13–${upper} qubits, only partitions observed in the Liège source library are listed. Reported numerical orders appear separately from rigorous bounds; no claim of complete coverage or proved numerical optimality is made. <a href="#sources">Library and provenance ↘</a></p>`;
  }
  function libraryRefinementSummary(refinement) {
    if(!refinement)return '';
    const number=scientific;
    return `<div class="library-refinement"><h4>Refined numerical vectors</h4><p class="source-description">${escape(refinement.primaryStates)} catalogue vectors${refinement.alternateStates?` and ${escape(refinement.alternateStates)} additional Dicke vectors`:''} were refined at their reported AC orders. The Frobenius-norm target is ${number(refinement.targetHSResidual)}; the largest residual after refinement is ${number(refinement.worstAfterHSResidual)}. AC verification remains numerical.</p><p class="micro-copy">Full-partition certificates apply to the saved catalogue vectors. Additional Dicke vectors are provided as downloads without a partition claim here.</p><div class="refinement-links">${refinement.downloadPath?`<a class="text-link" href="${escape(safeUrl(refinement.downloadPath))}">Download all refined coefficients ↗</a>`:''}${refinement.methodsPath?`<a class="text-link" href="${escape(safeUrl(refinement.methodsPath))}">Refinement methods ↗</a>`:''}${refinement.auditPath?`<a class="text-link" href="${escape(safeUrl(refinement.auditPath))}">Independent verification ↗</a>`:''}</div></div>`;
  }
  function renderSources(){
    const library=data.library;
    const libraryHTML=library?`<article class="source-item library-source"><span class="source-index">L</span><div><h3>Liège numerical libraries</h3><p class="source-description">${escape(library.provenanceNote||'Published filenames and links were recovered from saved copies of the primary pages. Numerical data were imported from matching local counterparts; current online byte identity is unverified because the live server timed out.')}</p>${library.validationNote?`<p class="source-description">${escape(library.validationNote)}</p>`:''}${libraryRefinementSummary(library.refinement)}<ul class="detail-sources">${(library.sourcePages||[]).map(page=>`<li><a href="${escape(safeUrl(page.url))}">${escape(page.title||page.url)} ↗</a></li>`).join('')}</ul>${library.missingFiles?.length?`<p class="micro-copy">Published files not available locally: ${library.missingFiles.map(file=>escape(typeof file==='string'?file:file.filename||file.id||'unavailable file')).join(', ')}.</p>`:''}</div></article>`:'';
    $('source-list').innerHTML=libraryHTML+(data.sources||[]).map((source,i)=>`<article class="source-item"><span class="source-index">${String(i+1).padStart(2,'0')}</span><div><h3><a href="${escape(sourceUrl(source))}">${escape(source.title||source.filename)} ↗</a></h3><p>${escape(Array.isArray(source.authors)?source.authors.join(', '):source.authors||'')}${source.year?' · '+escape(source.year):''}</p>${source.notes?`<p class="source-description">${escape(source.notes)}</p>`:''}<p class="source-filename">${escape(source.filename||'')}</p></div></article>`).join('');
    $('compiled-date').textContent=data.compiledDate?'Compiled '+data.compiledDate:'';
  }
  function revealProof(){if(location.hash.startsWith('#proof-')){const element=document.getElementById(location.hash.slice(1));if(element?.tagName==='DETAILS')element.open=true;}}
  $('view-majorana').addEventListener('click',openMajoranaViewer);
  $('close-majorana').addEventListener('click',()=>$('majorana-viewer').close());
  $('majorana-viewer').addEventListener('close',()=>{disposeViewerSphere();disposeViewerSphere=()=>{};document.body.classList.remove('majorana-viewer-open');const target=viewerReturnFocus?.isConnected&&!viewerReturnFocus.disabled?viewerReturnFocus:!$('view-majorana').disabled?$('view-majorana'):$('n-filter');target.focus({preventScroll:true});});
  $('majorana-viewer').addEventListener('click',event=>{if(event.target!==$('majorana-viewer'))return;const box=event.target.getBoundingClientRect();if(event.clientX<box.left||event.clientX>box.right||event.clientY<box.top||event.clientY>box.bottom)event.target.close();});
  $('viewer-witness').addEventListener('change',event=>{selected.stateId=event.target.value;renderDetail();updateHash();});
  $('state-table').addEventListener('click',event=>{const row=event.target.closest('tr[data-partition]');if(row)select(row.dataset.n,row.dataset.partition,null,true);});
  $('state-table').addEventListener('focusin',event=>{const button=event.target.closest('button[data-partition]');if(button)$('state-table').querySelectorAll('button[data-partition]').forEach(b=>b.tabIndex=b===button?0:-1);});
  $('state-table').addEventListener('keydown',event=>{
    const button=event.target.closest('button[data-partition]');if(!button)return;
    if(!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','Home','End','PageUp','PageDown'].includes(event.key))return;
    event.preventDefault();const buttons=[...$('state-table').querySelectorAll('button[data-partition]')];let index=buttons.indexOf(button);
    if(event.key==='Home')index=0;else if(event.key==='End')index=buttons.length-1;
    else index+=event.key==='PageUp'?-10:event.key==='PageDown'?10:['ArrowUp','ArrowLeft'].includes(event.key)?-1:1;
    buttons[Math.max(0,Math.min(buttons.length-1,index))]?.focus();
  });
  function refreshFilters(){ensureSelection();renderTable();renderDetail();updateHash();}
  $('status-filter').addEventListener('change',event=>{filter=event.target.value;refreshFilters();});
  $('n-filter').addEventListener('change',event=>{nFilter=event.target.value;$('state-table').scrollTop=0;refreshFilters();});
  document.addEventListener('click',event=>{const link=event.target.closest('a[href^="#proof-"]');if(link){const proof=document.getElementById(link.hash.slice(1));if(proof)proof.open=true;}});
  window.addEventListener('hashchange',()=>{if(location.hash.startsWith('#n=')){readHash();select(selected.n,selected.partitionKey,selected.stateId);}else revealProof();});
  readHash();ensureSelection();syncFilters();renderTable();renderDetail();renderProofs();renderLibraryScope();renderSources();typeset();revealProof();
  if(!location.hash||location.hash.startsWith('#n='))updateHash();
})();
