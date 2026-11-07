/* Offline AC-subspace catalogue. Mathematical data live in data/catalogue.js. */
(() => {
  'use strict';
  const data=window.AC_SUBSPACES;
  const $=id=>document.getElementById(id);
  const escape=value=>String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const math=(tex,display=false)=>`<span class="math" data-tex="${escape(tex)}" data-display="${display}">${escape(tex)}</span>`;
  function typeset(root=document){root.querySelectorAll('[data-tex]').forEach(el=>{if(window.katex)el.innerHTML=window.katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==='true',throwOnError:false,trust:false,strict:'ignore'});});}
  function safeUrl(value){const url=String(value||'').trim();return /^https?:\/\//i.test(url)||url&&!/^[a-z][a-z\d+.-]*:/i.test(url)&&!url.startsWith('//')?url:'#sources';}
  function references(refs=[]){return refs.map(ref=>`<li><a href="${escape(safeUrl(ref.url||ref.localUrl))}"${/^https?:/i.test(ref.url||'')?' target="_blank" rel="noopener"':''}>${escape(ref.label||ref.title||'Source')} ↗</a></li>`).join('');}
  function cleanBody(html=''){
    const template=document.createElement('template');template.innerHTML=html;
    const allowed=new Set(['P','UL','OL','LI','EM','STRONG','CODE','PRE','SPAN','DIV','A','SUB','SUP','BR','DETAILS','SUMMARY','BLOCKQUOTE','H3','H4','H5','TABLE','THEAD','TBODY','TR','TH','TD','CAPTION']);
    [...template.content.querySelectorAll('*')].forEach(el=>{if(!allowed.has(el.tagName)){el.replaceWith(document.createTextNode(el.textContent));return;}[...el.attributes].forEach(attr=>{if(!['class','data-tex','data-display','href','colspan','rowspan','scope'].includes(attr.name))el.removeAttribute(attr.name);if(attr.name==='href')el.setAttribute('href',safeUrl(attr.value));});});return template.innerHTML;
  }
  typeset();
  if(!data||!Array.isArray(data.entries)){$('subspace-table').innerHTML='<div class="empty-state"><h3>Source records are not loaded</h3><p>Keep the data folder beside this page to display the dimension catalogue.</p></div>';$('result-count').textContent='Source data unavailable';return;}
  const numbers=data.numbers||[...new Set(data.entries.map(e=>e.n))].sort((a,b)=>a-b);
  const orders=data.orders||[...new Set(data.entries.map(e=>e.t))].sort((a,b)=>a-b);
  const entries=new Map(data.entries.map(e=>[`${e.n}:${e.t}`,e]));
  const subspaces=new Map((data.subspaces||[]).map(s=>[s.id,s]));
  const proofs=new Map((data.proofs||[]).map(p=>[p.id,p]));
  const labels={proven:'Exact maximum',range:'Rigorous bounds',impossible:'No AC state'};
  let selected={n:7,t:2,subspaceId:null},filter='all',nFilter='all',toastTimer;
  const getEntry=()=>entries.get(`${selected.n}:${selected.t}`);
  const stateStatus=entry=>entry.upper===0?'impossible':Object.hasOwn(labels,entry.status)?entry.status:'range';
  const hasConjecture=entry=>Number.isInteger(entry.conjecturedDimension);
  const conjectureStatus=entry=>['confirmed','attained','open','superseded'].includes(entry.conjectureStatus)?entry.conjectureStatus:entry.lower>entry.conjecturedDimension?'superseded':entry.lower===entry.upper&&entry.lower===entry.conjecturedDimension?'confirmed':entry.lower>=entry.conjecturedDimension?'attained':'open';
  const conjectureLabels={confirmed:'Maximum confirmed',attained:'Dimension attained; maximum open',open:'Exact attainment still open',superseded:'Exceeded by an exact construction'};
  function matches(entry){
    if(nFilter!=='all'&&entry.n!==Number(nFilter))return false;
    if(filter==='all')return true;
    if(filter==='numerical')return Number.isInteger(entry.reportedDimension);
    if(filter==='conjectures')return hasConjecture(entry);
    if(filter==='conjectures-open')return hasConjecture(entry)&&['open','attained'].includes(conjectureStatus(entry));
    return stateStatus(entry)===filter;
  }
  const visible=()=>data.entries.filter(matches);
  function rangeText(entry){return entry.lower===entry.upper?String(entry.lower):`${entry.lower??0}–${entry.upper??'?'}`;}
  function rangeLatex(entry){return entry.lower===entry.upper?`k_{\\max}=${entry.lower}`:`${entry.lower??0}\\le k_{\\max}\\le ${entry.upper??'?'}`;}
  function hashFor(){const params=new URLSearchParams({n:selected.n,t:selected.t});if(selected.subspaceId)params.set('subspace',selected.subspaceId);return '#'+params.toString();}
  function updateHash(){const hash=hashFor();if(location.hash!==hash)history.replaceState(null,'',hash);}
  function readHash(){const params=new URLSearchParams(location.hash.slice(1)),n=Number(params.get('n')),t=Number(params.get('t'));if(entries.has(`${n}:${t}`))selected={n,t,subspaceId:params.get('subspace')};}
  function ensureSelection(){if(!getEntry()||!matches(getEntry())){const first=visible()[0];if(first)selected={n:first.n,t:first.t,subspaceId:null};}}
  function select(n,t,subspaceId=null,focus=false){selected={n:Number(n),t:Number(t),subspaceId};if(!entries.has(`${selected.n}:${selected.t}`))return;if(!matches(getEntry())){filter='all';nFilter='all';}$('status-filter').value=filter;$('n-filter').value=nFilter;renderMatrix();renderDetail();updateHash();if(focus)$('subspace-table').querySelector(`[data-n="${selected.n}"][data-t="${selected.t}"]`)?.focus({preventScroll:true});}
  function notify(message){$('toast').textContent=message;$('toast').classList.add('visible');clearTimeout(toastTimer);toastTimer=setTimeout(()=>$('toast').classList.remove('visible'),2600);}
  function renderMatrix(){
    const ns=nFilter==='all'?numbers:numbers.filter(n=>n===Number(nFilter));let count=0;
    const rows=ns.map(n=>`<tr><th scope="row">${n}</th>${orders.map(t=>{
      const entry=entries.get(`${n}:${t}`);if(!entry)return '<td class="outside-domain" aria-label="No catalogue record">—</td>';
      if(!matches(entry))return '<td class="filtered-cell" aria-label="Hidden by evidence filter">·</td>';
      count++;const status=stateStatus(entry),isSelected=n===selected.n&&t===selected.t,unresolved=status==='range',conjectured=unresolved&&hasConjecture(entry),reported=unresolved&&Number.isInteger(entry.reportedDimension)&&(!conjectured||entry.reportedDimension!==entry.conjecturedDimension);
      const accessible=`${n} ${n===1?'qubit':'qubits'}, order ${t}: dimension ${rangeText(entry)}, ${labels[status]}${conjectured?`, paper conjectures maximum ${entry.conjecturedDimension}: ${conjectureLabels[conjectureStatus(entry)]}`:''}${reported?`, numerical achieved dimension ${entry.reportedDimension}`:''}`;
      return `<td><button type="button" class="state-cell ${status}${isSelected?' selected':''}" data-n="${n}" data-t="${t}" aria-pressed="${isSelected}" aria-label="${escape(accessible)}"><span class="dimension-value">${escape(rangeText(entry))}</span><span class="cell-label">${status==='proven'?'exact maximum':status==='impossible'?'no AC state':'rigorous bounds'}</span>${conjectured?`<span class="conjecture-cell" title="Paper’s conjectured maximum: ${entry.conjecturedDimension} · ${escape(conjectureLabels[conjectureStatus(entry)])}">Conj. ${entry.conjecturedDimension}</span>`:''}${reported?`<span class="numeric-cell" title="Numerically achieved dimension">Num. ${entry.reportedDimension}</span>`:''}</button></td>`;
    }).join('')}</tr>`).join('');
    $('subspace-table').innerHTML=`<table class="state-grid"><caption class="sr-only">Largest dimension of anticoherent subspaces by qubit number n and order t. Paper conjectures and numerical achieved dimensions are separate from rigorous bounds.</caption><thead><tr><th scope="col">${math('n\\backslash t')}</th>${orders.map(t=>`<th scope="col">${t}</th>`).join('')}</tr></thead><tbody>${rows}</tbody></table>${count?'':'<p class="no-matches">No cells match these filters.</p>'}`;
    $('result-count').textContent=`${count} ${count===1?'cell':'cells'} · ${subspaces.size} explicit ${subspaces.size===1?'basis':'bases'}`;typeset($('subspace-table'));
  }
  function decimal(value){const [re,im]=Array.isArray(value)?value:[Number(value),0];const f=x=>Math.abs(x)<1e-14?'0':Number(x).toPrecision(9).replace(/(?:\.0+|(\.\d*?)0+)$/,'$1');return Math.abs(im)<1e-14?f(re):Math.abs(re)<1e-14?`${f(im)}i`:`${f(re)} ${im<0?'−':'+'} ${f(Math.abs(im))}i`;}
  function scientific(value){if(value==null)return 'Not supplied';const n=Number(value);if(!Number.isFinite(n))return escape(value);if(n===0)return '0';if(n>=.001&&n<1000)return escape(Number(n.toPrecision(5)));const [a,b]=n.toExponential(3).split('e');return math(`${a}\\times10^{${Number(b)}}`);}
  // A small display parser for the saved exact expressions. It never evaluates code.
  function exactLatex(value){
    const source=String(value??''),tokens=[];let offset=0;
    while(offset<source.length){if(/^\s*$/.test(source.slice(offset)))break;const match=/^\s*(sqrt|I|\d+(?:\.\d+)?|\*\*|[+\-*/()])/.exec(source.slice(offset));if(!match)return null;tokens.push(match[1]);offset+=match[0].length;}
    let i=0;
    function primary(){const token=tokens[i++];if(token==='sqrt'){if(tokens[i++]!=='(')throw Error();const inside=addition();if(tokens[i++]!==')')throw Error();return {type:'sqrt',inside};}if(token==='('){const inside=addition();if(tokens[i++]!==')')throw Error();return inside;}if(token==='I')return {type:'atom',text:'i'};if(token&&/^\d/.test(token))return {type:'atom',text:token};throw Error();}
    function power(){let left=primary();if(tokens[i]==='**'){i++;left={type:'pow',left,right:unary()};}return left;}
    function unary(){if(tokens[i]==='+'||tokens[i]==='-'){const op=tokens[i++];return {type:'unary',op,inside:unary()};}return power();}
    function product(){let left=unary();while(tokens[i]==='*'||tokens[i]==='/'){const op=tokens[i++];left={type:op==='*'?'mul':'div',left,right:unary()};}return left;}
    function addition(){let left=product();while(tokens[i]==='+'||tokens[i]==='-'){const op=tokens[i++];left={type:'add',op,left,right:product()};}return left;}
    function render(node,precedence=0){let tex,p;if(node.type==='atom'){tex=node.text;p=5;}else if(node.type==='sqrt'){tex='\\sqrt{'+render(node.inside)+'}';p=5;}else if(node.type==='div'){tex='\\frac{'+render(node.left)+'}{'+render(node.right)+'}';p=5;}else if(node.type==='pow'){tex='{'+render(node.left,4)+'}^{'+render(node.right)+'}';p=4;}else if(node.type==='unary'){tex=node.op+render(node.inside,3);p=3;}else if(node.type==='mul'){tex=render(node.left,2)+'\\,'+render(node.right,2);p=2;}else{tex=render(node.left,1)+node.op+render(node.right,node.op==='-'?2:1);p=1;}return p<precedence?'\\left('+tex+'\\right)':tex;}
    try{const tree=addition();return i===tokens.length?render(tree):null;}catch{return null;}
  }
  function derivationDetails(subspace){
    const d=subspace.derivation;if(!d)return '';
    const vector=(items)=>'\\left('+items.map(x=>exactLatex(x)||String(x)).join(',\\,')+'\\right)';
    let html='<details class="coefficient-disclosure derivation-details"><summary>Exact reconstruction: weights & phases</summary>';
    if(d.argument)html+=`<p class="detail-note">${escape(d.argument)}</p>`;
    (d.supports||[]).forEach((support,index)=>{html+=`<div class="derivation-equation">${math('R_{'+(index+1)+'}='+vector(support),true)}</div>`;if(d.probabilities?.[index])html+=`<div class="derivation-equation">${math('p_{'+(index+1)+'}='+vector(d.probabilities[index]),true)}</div>`;});
    if(d.probabilities?.length)html+='<p class="micro-copy">Each probability list follows the corresponding excitation support. Components are the square roots of these probabilities multiplied by the phases below.</p>';
    const phaseRows=d.phaseChoice||d.phases||(d.A_phases?{A:d.A_phases,B:d.B_phases}:{});
    for(const [name,phases] of Object.entries(phaseRows)){if(Array.isArray(phases))html+=`<div class="derivation-equation">${math('\\boldsymbol{\\eta}_{'+name+'}='+vector(phases),true)}</div>`;}
    const cosine=d.phaseDifferenceCosine||d.phaseCosine;
    if(cosine){const tex=exactLatex(cosine);if(tex)html+=`<div class="derivation-equation">${math('\\cos\\phi='+tex,true)}</div>`;}
    if(d.commonPhase){const tex=exactLatex(d.commonPhase);if(tex)html+=`<div class="derivation-equation">${math('e^{i\\phi}='+tex,true)}</div>`;}
    return html+'</details>';
  }
  function coefficientTable(vector){
    const implicit=!!vector.implicitAlgebraicDefinition;
    const count=Math.max(vector.coefficients?.length||0,vector.exactCoefficients?.length||0);
    if(!count)return '<p class="micro-copy">Component data are not attached to this basis vector.</p>';
    const support=vector.support||Array.from({length:count},(_,r)=>r).filter(r=>vector.exactCoefficients?.[r]!=null?String(vector.exactCoefficients[r])!=='0':Math.hypot(...(Array.isArray(vector.coefficients?.[r])?vector.coefficients[r]:[Number(vector.coefficients?.[r]||0),0]))>1e-14);
    return `<details class="coefficient-disclosure"><summary>Dicke coefficients <span>${support.length} nonzero</span></summary><div class="coefficient-scroll"><table class="coefficient-table"><thead><tr><th scope="col">r</th>${implicit?'':'<th scope="col">Exact amplitude</th>'}<th scope="col">${implicit?'Approximate amplitude':'Decimal'}</th></tr></thead><tbody>${support.map(r=>`<tr><th scope="row">${r}</th>${implicit?'':`<td class="exact-amplitude">${vector.coefficientLatex?.[r]!=null?math(vector.coefficientLatex[r]):vector.exactCoefficients?.[r]!=null?(exactLatex(vector.exactCoefficients[r])!=null?math(exactLatex(vector.exactCoefficients[r])):escape(vector.exactCoefficients[r])):'—'}</td>`}<td class="decimal-amplitude">${vector.coefficients?.[r]!=null?escape(decimal(vector.coefficients[r])):'—'}</td></tr>`).join('')}</tbody></table></div><p class="micro-copy">Normalized Dicke basis; r counts excitations. Decimal components are approximations.${implicit?' The linked certificate defines the exact algebraic state.':''}</p></details>`;
  }
  function majoranaRecord(vector,subspace,index){return {id:vector.id||`${subspace.id}-basis-${index+1}`,n:vector.n||subspace.n,label:`Basis vector ${index+1} · ${vector.label||subspace.label||subspace.id}`,points:vector.points||vector.majoranaPoints||[],evidence:`Majorana representation of this basis vector only; it does not describe the entire subspace. ${vector.evidence||'Coordinates are a numerical rendering of the saved Dicke components.'}`};}
  // Resize coefficient parentheses without changing Dicke indices or existing sized delimiters.
  function stateFormula(tex){return String(tex).replace(/(^|[=+\-])\(([\s\S]*?)\)(?=\|D_)/g,'$1\\left($2\\right)');}
  function basisVector(vector,subspace,index){
    const record=majoranaRecord(vector,subspace,index),hasPoints=record.points.length>0;
    const group=hasPoints?window.MajoranaViewer?.groupSummary?.(record,true)||'':'';
    return `<details class="basis-vector"${index===0?' open':''}><summary><span class="basis-index">${math(`|u_{${index+1}}\\rangle`)}</span><span>${escape(vector.label||`Basis vector ${index+1}`)}</span></summary>${vector.stateLatex?`<div class="state-formula" role="region" tabindex="0" aria-label="Dicke formula for basis vector ${index+1}">${math(stateFormula(vector.stateLatex),true)}</div>`:''}${vector.coefficientDefinitionLatex?`<div class="state-formula" role="region" tabindex="0" aria-label="Coefficient definitions for basis vector ${index+1}">${math(vector.coefficientDefinitionLatex,true)}</div>`:''}${group}${hasPoints?`<button type="button" class="view-majorana-button" data-basis="${index}" aria-haspopup="dialog" aria-controls="majorana-viewer">View this vector’s Majorana sphere ↗</button><p class="basis-point-note">Point groups and symmetry axes refer to this vector, not the subspace.</p>`:'<p class="basis-point-note">Majorana coordinates are not attached to this vector.</p>'}${coefficientTable(vector)}</details>`;
  }
  function conjectureCard(entry){
    if(!hasConjecture(entry))return '';
    const status=conjectureStatus(entry),value=entry.conjecturedDimension;
    const meaning=status==='confirmed'?(value===0?'The rigorous bounds confirm that no nonzero AC subspace exists.':'Matching rigorous lower and upper bounds now confirm the paper’s conjectured maximum.'):status==='attained'?(value===0?'The conjectured exclusion remains open; no matching zero upper bound is established.':'An exact construction attains this dimension. Proving that no larger subspace exists remains open.'):status==='superseded'?'An exact construction exceeds this proposed value, so it is no longer a possible maximum.':'The rigorous lower bound is below this proposed value. Exact attainment and maximality have not been established here.';
    return `<div class="conjectured-dimension ${status}"><h4>Paper’s conjectured maximum</h4><div class="conjecture-value">${math('k_{\\max}='+value)}</div><span class="conjecture-status">${escape(conjectureLabels[status])}</span><p>${escape(meaning)}</p>${entry.conjectureReferences?.length?`<ul class="detail-sources">${references(entry.conjectureReferences)}</ul>`:''}</div>`;
  }
  function renderDetail(){
    window.MajoranaViewer?.close();window.SubspaceViewer?.close();const entry=getEntry();
    if(!entry||!matches(entry)){$('result-detail').innerHTML='<p class="eyebrow">Selected pair</p><p>No cells match the current filters.</p>';return;}
    const status=stateStatus(entry),candidates=(entry.subspaceIds||[]).map(id=>subspaces.get(id)).filter(s=>s&&s.n===entry.n&&s.t>=entry.t),subspace=candidates.find(s=>s.id===selected.subspaceId)||candidates[0];selected.subspaceId=subspace?.id||null;
    let html=`<div class="detail-heading"><p class="eyebrow">Selected pair</p><span class="badge ${status}">${labels[status]}</span></div><h3>${entry.n} ${entry.n===1?'qubit':'qubits'} <span class="detail-separator">/</span> order ${entry.t}</h3><div class="dimension-metric">${math(rangeLatex(entry),true)}</div><p class="dimension-description">${status==='impossible'?'No nonzero subspace satisfies this AC order.':status==='proven'?'The lower and upper dimension bounds coincide.':entry.lower>0?'An exact construction or theorem supplies the lower bound; a general exclusion supplies the upper bound.':'No positive rigorous lower bound is recorded; the upper bound limits every possible AC subspace.'}</p>`;
    html+=conjectureCard(entry);
    if(Number.isInteger(entry.reportedDimension))html+=`<div class="reported-dimension"><strong>Numerical example: dimension ${entry.reportedDimension}</strong><p>A numerically achieved dimension, separate from the rigorous bounds and the paper’s proposed maximum. A small residual is not an exact certificate.</p>${entry.reportedReferences?.length?`<ul class="detail-sources">${references(entry.reportedReferences)}</ul>`:''}</div>`;
    if(entry.numericalSearch){
      const survey=entry.numericalSearch;
      const message=survey.noCandidateFound?'No candidate reached the acceptance threshold. This does not prove nonexistence.':survey.failedLargerDimensions.length?`No candidate was found at the tested larger ${survey.failedLargerDimensions.length===1?'dimension':'dimensions'} ${survey.failedLargerDimensions.join(', ')}. Larger subspaces remain unexcluded unless ruled out by the rigorous bound.`:survey.matchesRigorousUpper&&status!=='proven'?'The numerical dimension matches the rigorous upper bound. An exact existence certificate is still required to prove equality.':'The survey retains the rigorous bounds and records numerical evidence separately.';
      html+=`<div class="reported-dimension"><strong>Numerical survey · N = 21–24</strong><p>${escape(message)}</p><p class="micro-copy">${survey.attemptCount?`${survey.attemptCount} recorded search starts for this order. `:''}<a href="${escape(survey.reportUrl)}">Search results and residuals ↗</a></p></div>`;
    }
    if(entry.notes)html+=`<p class="detail-note">${escape(entry.notes)}</p>`;
    if(candidates.length){
      html+=`<label class="witness-label" for="subspace-select">Representative subspace<select id="subspace-select">${candidates.map(s=>`<option value="${escape(s.id)}"${s.id===subspace.id?' selected':''}>${escape(s.label||s.id)} · k = ${s.dimension} · ${s.evidence==='exact'?'exact':'numerical'}</option>`).join('')}</select></label><div class="subspace-metadata"><span>Dimension ${subspace.dimension}</span><span>AC order ${subspace.t}</span><span>${subspace.evidence==='exact'?'Exact construction':'Numerical construction'}</span></div>`;
      if(subspace.notes)html+=`<p class="detail-note">${escape(subspace.notes)}</p>`;
      if(subspace.exactnessKind==='isolated-algebraic-root'){
        const proofUrl=subspace.references?.find(ref=>/\.html(?:#|$)/.test(ref.url))?.url||`research/algebraic-k3/index.html#n${subspace.n}`;
        html+=`<div class="reported-dimension"><strong>Exact algebraic coefficients</strong><p>The named parameters are the unique real solution of a rational polynomial system inside a certified isolating box. Existence and all whole-subspace identities are proved exactly; the decimals are approximations of that root.</p><a href="${escape(proofUrl)}">Defining equations, parameter box and proof ↗</a></div>`;
      }
      const verification=subspace.verification||{};
      html+=`<p class="subspace-evidence">${verification.exact===true?'The basis orthonormality and all cross-reduction identities are verified exactly.':subspace.evidence==='exact'?'This is an exact source construction; see the attached evidence.':'The basis and its cross reductions are checked numerically; the reported dimension is not a rigorous lower bound.'}</p>`;
      if(verification.method)html+=`<p class="micro-copy">${escape(verification.method)}</p>`;
      if(verification.orthogonalityError!=null||verification.crossReductionError!=null)html+=`${verification.exact===true?'<p class="micro-copy">Floating-point cross-check of the displayed decimal components; the exact identities above have zero residual.</p>':''}<dl class="subspace-checks">${verification.orthogonalityError!=null?`<div><dt>Orthogonality error</dt><dd>${scientific(verification.orthogonalityError)}</dd></div>`:''}${verification.crossReductionError!=null?`<div><dt>Cross-reduction error</dt><dd>${scientific(verification.crossReductionError)}</dd></div>`:''}</dl>`;
      html+=`<div class="whole-subspace-entry"><p class="eyebrow">Whole subspace</p><p>Explore the basis-independent projector: its Husimi function, rotational symmetry axes and operator Majorana constellations.</p><button type="button" id="view-subspace" class="view-majorana-button" aria-haspopup="dialog" aria-controls="subspace-viewer">View whole subspace ↗</button></div>`;
      html+=derivationDetails(subspace);
      html+=`<h4 class="basis-heading">${subspace.evidence==='exact'?'An orthonormal basis':'A numerical basis'}</h4>${subspace.evidence==='numerical'?'<p class="micro-copy">The displayed decimal vectors satisfy orthogonality and the cross-reduction identities only to the reported numerical accuracy.</p>':''}<p class="basis-warning">Each sphere shows one basis vector. A change of basis can change these constellations while preserving the same subspace.</p>`;
      if(subspace.derivationLatex)html+=`<div class="derivation-equation basis-parameter" aria-label="Basis parameter definition">${math(subspace.derivationLatex,true)}</div>`;
      if(subspace.parameterLatex)html+=`<div class="derivation-equation basis-parameter" aria-label="Algebraic parameter coordinates">${math(subspace.parameterLatex,true)}</div>`;
      html+=(subspace.basis||[]).map((vector,index)=>basisVector(vector,subspace,index)).join('')||'<p class="micro-copy">The basis components are not attached to this record.</p>';
      if(subspace.references?.length)html+=`<h4 class="detail-label">Construction evidence</h4><ul class="detail-sources">${references(subspace.references)}</ul>`;
    }else html+=`<p class="detail-note">${status==='impossible'?'There is no representative to visualize.':'No explicit basis is attached to this cell. The bounds above remain the recorded result.'}</p>`;
    const linkedProofs=(entry.proofIds||[]).map(id=>proofs.get(id)).filter(Boolean);
    if(linkedProofs.length)html+=`<h4 class="detail-label">Bounds & proof records</h4><ul class="detail-proof-links">${linkedProofs.map(p=>`<li><a href="#proof-${escape(p.id)}">${escape(p.title)} ↘</a></li>`).join('')}</ul>`;
    html+='<div class="detail-actions"><button type="button" id="download-record">Download record ↓</button><button type="button" id="copy-link">Copy link ↗</button></div>';
    $('result-detail').innerHTML=html;typeset($('result-detail'));
    $('view-subspace')?.addEventListener('click',()=>window.SubspaceViewer?.open(subspace));
    $('subspace-select')?.addEventListener('change',event=>{selected.subspaceId=event.target.value;renderDetail();updateHash();$('subspace-select')?.focus({preventScroll:true});});
    $('result-detail').querySelectorAll('button[data-basis]').forEach(button=>button.addEventListener('click',()=>{const index=Number(button.dataset.basis),vector=subspace.basis[index];window.MajoranaViewer?.open(majoranaRecord(vector,subspace,index));}));
    $('download-record').addEventListener('click',()=>download({entry,selectedSubspace:subspace||null,proofs:linkedProofs},`ac-subspace-n${entry.n}-t${entry.t}.json`));
    $('copy-link').addEventListener('click',async()=>{const url=location.href.split('#')[0]+hashFor();try{if(!navigator.clipboard?.writeText)throw new Error();await navigator.clipboard.writeText(url);notify('Link copied');}catch{const input=document.createElement('textarea');input.value=url;input.className='copy-fallback';document.body.append(input);input.select();const ok=document.execCommand('copy');input.remove();notify(ok?'Link copied':'Copy the address from your browser to share this result.');}});
  }
  function download(payload,filename){const url=URL.createObjectURL(new Blob([JSON.stringify(payload,null,2)+'\n'],{type:'application/json'})),link=document.createElement('a');link.href=url;link.download=filename;link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);notify('Record downloaded');}
  function renderProofs(){
    $('proof-list').innerHTML=(data.proofs||[]).map(proof=>`<details class="proof-item" id="proof-${escape(proof.id)}"><summary><span>${escape(proof.title)}</span><span class="proof-expand" aria-hidden="true">+</span></summary><div class="proof-body">${proof.summary?`<p>${escape(proof.summary)}</p>`:''}${proof.formulaLatex?math(proof.formulaLatex,true):''}${cleanBody(proof.bodyHtml||'')}${proof.references?.length?`<ul class="detail-sources">${references(proof.references)}</ul>`:''}</div></details>`).join('')||'<p class="section-intro">No proof records are attached.</p>';
  }
  function renderSources(){
    const caveats=data.sourceCaveats||[];$('source-caveats').innerHTML=caveats.length?`<details class="caveat-box"><summary>Scope & reading notes <span>${caveats.length}</span></summary><div>${caveats.map(note=>`<p>${escape(typeof note==='string'?note:note.notes||note.issue||'')}</p>`).join('')}</div></details>`:'';
    $('source-list').innerHTML=(data.sources||[]).map((source,index)=>`<article class="source-item" id="source-${escape(source.id||index)}"><span class="source-index">${String(index+1).padStart(2,'0')}</span><div><h3>${escape(source.title||source.id)}</h3><p>${escape(Array.isArray(source.authors)?source.authors.join(', '):source.authors||'')}${source.year?' · '+escape(source.year):''}</p>${source.notes?`<p class="source-description">${escape(source.notes)}</p>`:''}<p class="source-actions">${source.localUrl?`<a href="${escape(safeUrl(source.localUrl))}">Local source ↗</a>`:''}${source.url&&source.url!==source.localUrl?`<a href="${escape(safeUrl(source.url))}" target="_blank" rel="noopener">Source website ↗</a>`:''}</p></div></article>`).join('');
    $('coverage-note').textContent=typeof data.scope==='string'?data.scope:'';
  }
  function revealProof(){if(location.hash.startsWith('#proof-')){const target=document.getElementById(location.hash.slice(1));if(target?.tagName==='DETAILS')target.open=true;}}
  $('n-filter').innerHTML='<option value="all">All</option>'+numbers.map(n=>`<option value="${n}">${n}</option>`).join('');
  $('n-filter').addEventListener('change',event=>{nFilter=event.target.value;ensureSelection();renderMatrix();renderDetail();updateHash();});
  $('status-filter').addEventListener('change',event=>{filter=event.target.value;ensureSelection();renderMatrix();renderDetail();updateHash();});
  $('subspace-table').addEventListener('click',event=>{const button=event.target.closest('button[data-n]');if(button)select(button.dataset.n,button.dataset.t,null,true);});
  $('subspace-table').addEventListener('keydown',event=>{const button=event.target.closest('button[data-n]');if(!button||!['ArrowLeft','ArrowRight','ArrowUp','ArrowDown','Home','End'].includes(event.key))return;event.preventDefault();const ns=nFilter==='all'?numbers:numbers.filter(n=>n===Number(nFilter));let row=ns.indexOf(Number(button.dataset.n)),col=orders.indexOf(Number(button.dataset.t));if(event.key==='Home')col=0;else if(event.key==='End')col=orders.length-1;else{row+=event.key==='ArrowDown'?1:event.key==='ArrowUp'?-1:0;col+=event.key==='ArrowRight'?1:event.key==='ArrowLeft'?-1:0;}if(row<0||row>=ns.length||col<0||col>=orders.length)return;$('subspace-table').querySelector(`[data-n="${ns[row]}"][data-t="${orders[col]}"]`)?.focus();});
  document.addEventListener('click',event=>{const link=event.target.closest('a[href^="#proof-"]');if(link){const target=document.getElementById(link.hash.slice(1));if(target)target.open=true;}});
  window.addEventListener('hashchange',()=>{if(location.hash.startsWith('#n=')){readHash();select(selected.n,selected.t,selected.subspaceId);}else revealProof();});
  readHash();ensureSelection();renderMatrix();renderDetail();renderProofs();renderSources();typeset();revealProof();if(!location.hash||location.hash.startsWith('#n='))updateHash();
})();
