(function(){
  'use strict';
  const $=id=>document.getElementById(id), esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const tex=(source)=>window.katex?katex.renderToString(source,{throwOnError:false,strict:false}):esc(source);
  const names={2:'Two',3:'Three',4:'Four',5:'Five'};
  let n=3,data,named,selected,current,sector,sphere,previewTimer,activeCurve,interactiveState=null,previewing=false,pointDragging=false;
  const fmt=(v,d=6)=>Math.abs(v)<1e-10?'0':Math.abs(v-Math.round(v))<1e-9?String(Math.round(v)):Number(v.toFixed(d)).toString();
  const title=s=>s.label||((s.family||'Pure symmetric')+' state');
  function complex(c){const [r,i]=c;return Math.abs(i)<1e-10?fmt(r):Math.abs(r)<1e-10?(Math.abs(i-1)<1e-9?'i':Math.abs(i+1)<1e-9?'−i':fmt(i)+'i'):`${fmt(r)} ${i<0?'−':'+'} ${fmt(Math.abs(i))}i`;}
  function formula(c){const terms=[];c.forEach((a,k)=>{if(Math.hypot(...a)<1e-9)return;let x=complex(a).replace(/−/g,'-');const coefficient=x==='1'?'':x==='-1'?'-':(Math.abs(a[0])>1e-9&&Math.abs(a[1])>1e-9?'('+x+')':x);terms.push(coefficient+`\\,|D_{${n}}^{(${k})}\\rangle`);});return terms.map((t,i)=>'<span class="term">'+tex((i&&t[0]!=='-'?'+ ':'')+t)+'</span>').join('');}
  function parameterBlocks(source){
    const at=source.indexOf('\\begin{aligned}');
    if(at<0)return source.split(/,?\\qquad\s*|,\\quad\s*/);
    const prefix=source.slice(0,at).replace(/(?:,?\s*\\(?:qquad|quad)\s*)+$/,'').trim();
    return prefix?[prefix,source.slice(at)]:[source];
  }
  function paintMath(el,source){if(el._mathSource===source&&el.innerHTML===el._mathHTML)return;el._mathSource=source;el.innerHTML=AtlasMath.html(source);el._mathHTML=el.innerHTML;}
  function showExpression(s,coefficients,{live=false}={}){const exact=s.analyticalExpression||(window.EXACT_SMALL_EXPRESSIONS||{})[s.id];current.analyticalExpression=exact||null;if(exact)paintMath($('state-formula'),'$\\displaystyle '+exact.latex+'$');else{$('state-formula')._mathSource=null;$('state-formula').innerHTML=formula(coefficients);}const inline=!!exact?.parameterLatex&&!exact.parameterLatex.includes('begin{aligned}')&&exact.parameterLatex.length<130;$('parameter-inline').hidden=!inline;$('parameter-inline').innerHTML=inline?AtlasMath.html('$'+exact.parameterLatex+'$'):'';const details=$('analytical-details'),stateKey=s.id||'unnamed';if(details.dataset.stateId!==stateKey){details.open=false;details.dataset.stateId=stateKey;}details.hidden=!exact?.parameterLatex||inline;if(exact?.parameterLatex&&!inline&&!live){const source=exact.parameterLatex;const parts=parameterBlocks(source);$('analytical-parameters').innerHTML=parts.map(t=>'<div class="parameter-equation">'+AtlasMath.html('$'+t+'$')+'</div>').join('');}$('dicke-shorthand').hidden=true;}
  function fitSectorValues(){
    const grid=$('sector-values');grid.classList.remove('exact-2','exact-1');
    const entries=[...grid.querySelectorAll('.sector-value')],width=grid.clientWidth;
    const need=Math.max(0,...entries.map(el=>(el.querySelector('.sector-symbol')?.getBoundingClientRect().width||0)+(el.querySelector('.sector-number>.katex')?.getBoundingClientRect().width||0)+12));
    if(need>(width-14)/3)grid.classList.add(need>(width-7)/2?'exact-1':'exact-2');
  }
  function renderPurityValues(state){
    const {purities}=state;
    const exact=PurityCoordinates.exactPurityLabels(n,state.exactSectors);
    $('sector-values').innerHTML=purities.map((v,k)=>{
      const relation=exact[k]?'=':'\\approx';
      return '<div class="sector-value"><span class="sector-symbol">'+tex('P_{'+(k+1)+'}'+relation)+'</span><span class="sector-number">'+tex(exact[k]??fmt(v))+'</span></div>';
    }).join('');
    requestAnimationFrame(fitSectorValues);
  }
  function stateData(s){
    const coefficients=SectorMath.canonicalState(s.coefficients);
    const sectors=SectorMath.sectorLengths(coefficients);
    const purities=Array.from({length:n},(_,k)=>k===n-1?1:SectorMath.reducedPurity(coefficients,k+1));
    return {...s,N:n,coefficients,sectors,purities};
  }
  function intervalState(x,size=n){
    if(!Number.isFinite(x))return null;
    const s=SectorConstructors[size===2?"n2":"n3"](Number(x.toPrecision(12)));
    if(!s)return null;
    [x]=s.sectors;
    return {...s,id:`n${size}-point:${x}`,interactivePoint:true,
      label:names[size]+'-qubit interval state',mathLabel:'$|\\psi(x)\\rangle$',
      analyticalExpression:{
        latex:'\\sqrt{\\frac{1+r}{2}}\\,|D_{'+size+'}^{(0)}\\rangle+\\sqrt{\\frac{1-r}{2}}\\,|D_{'+size+'}^{('+size+')}\\rangle',
        parameterLatex:`x=${x},\\qquad r=\\sqrt{x/${size}}`,
        analyticalStatus:'Exact generalized GHZ family covering the entire interval; clicked coordinate stored at 12 significant digits'
      }};
  }
  function intervalStateFromId(id){
    const m=/^n([23])-point:([^:]+)$/.exec(id||'');
    return m&&Number(m[1])===n?intervalState(Number(m[2])):null;
  }
  function triangleState(x,y){
    if(!Number.isFinite(x)||!Number.isFinite(y))return null;
    const s=SectorConstructors.n4(Number(x.toPrecision(12)),Number(y.toPrecision(12)));
    if(!s)return null;
    [x,y]=s.sectors;
    return {...s,id:`n4-point:${x}:${y}`,interactivePoint:true,
      label:'Triangle state',mathLabel:'$|\\psi(x,y)\\rangle$',
      analyticalExpression:{
        latex:'\\cos\\frac\\theta2\\left(\\cos\\frac\\phi2\\,|D_4^{(0)}\\rangle+\\sin\\frac\\phi2\\,|D_4^{(4)}\\rangle\\right)+i\\sin\\frac\\theta2\\,|D_4^{(2)}\\rangle',
        parameterLatex:`x=${x},\\quad y=${y},\\qquad\\begin{aligned}Y&=3-y/2,\\quad t=\\sqrt{4-x-2Y},\\\\b&=\\sqrt{\\frac{Y}{2(2+t)}},\\quad a=\\sqrt{1-b^2},\\\\\\theta&=2\\arctan(b/a),\\\\\\phi&=\\arcsin\\!\\left(1-\\frac{x}{4+2t-Y}\\right).\\end{aligned}`,
        analyticalStatus:'Exact entire-triangle inverse; clicked coordinates stored at 12 significant digits'
      }};
  }
  function triangleStateFromId(id){
    const m=/^n4-point:([^:]+):([^:]+)$/.exec(id||'');
    return m?triangleState(Number(m[1]),Number(m[2])):null;
  }
  function fiveQubitState(x,y){
    if(!Number.isFinite(x)||!Number.isFinite(y))return null;
    const s=SectorConstructors.n5(Number(x.toPrecision(12)),Number(y.toPrecision(12)));
    if(!s)return null;
    [x,y]=s.sectors;
    const p=s.parameter,xy=`x=${x},\\quad y=${y}`,value=v=>Number(v.toPrecision(12));
    let latex,definitions;
    if(p.branch==='parity'){
      latex='-\\sqrt{\\frac{R-t}{2}}|D_5^{(0)}\\rangle+\\sqrt t\\,|D_5^{(2)}\\rangle+\\sqrt{\\frac{2-R-t}{2}}|D_5^{(4)}\\rangle';
      definitions=`R=\\frac{3-\\sqrt{5x}}4,\\quad t\\approx${value(p.t)},\\qquad\\begin{aligned}q&=2R-2/3-4t,\\\\\\zeta&=-2\\sqrt{5t(R-t)}+6\\sqrt{t(2-R-t)},\\\\y&=\\frac{10}3+\\frac35q^2+\\frac15\\zeta^2.\\end{aligned}`;
    }else if(p.branch==='cat'){
      latex='\\sqrt{1-b-c}\\,|D_5^{(0)}\\rangle+\\sqrt b\\,|D_5^{(1)}\\rangle+\\sqrt c\\,|D_5^{(5)}\\rangle';
      definitions=`b\\approx${value(p.b)},\\qquad\\begin{aligned}c&=\\frac{5-b-\\sqrt{17b^2-10b+5x}}{10},\\\\y&=10-\\frac{32}5b^2-16bc.\\end{aligned}`;
    }else{
      latex='\\cos\\alpha\\,|\\chi\\rangle+e^{i\\delta}\\sin\\alpha\\,\\Theta|\\chi\\rangle';
      definitions=`t\\approx${value(p.t)},\\quad\\epsilon=${p.sign},\\qquad\\begin{aligned}|\\chi\\rangle&=a|D_5^{(0)}\\rangle+b|D_5^{(2)}\\rangle+c|D_5^{(4)}\\rangle,\\\\a&=\\epsilon\\sqrt{\\frac{3/4-t}{2}},\\quad b=\\sqrt t,\\quad c=\\sqrt{\\frac{5/4-t}{2}},\\\\\\Theta|\\chi\\rangle&=c|D_5^{(1)}\\rangle+b|D_5^{(3)}\\rangle+a|D_5^{(5)}\\rangle,\\\\\\alpha&=\\frac12\\arcsin\\frac{\\sqrt{5x}}{2\\kappa},\\quad\\kappa\\approx${value(p.kappa)},\\\\\\delta&=${p.axis==='x'?'0':'\\pi/2'}.\\end{aligned}`;
    }
    return {...s,id:`n5-point:${x}:${y}`,interactivePoint:true,label:'Five-qubit region state',mathLabel:'$|\\psi(x,y)\\rangle$',
      analyticalExpression:{latex,parameterLatex:xy+',\\qquad '+definitions,analyticalStatus:'Analytically bracketed exact state family; parameters evaluated numerically'}};
  }
  function fiveQubitStateFromId(id){
    const m=/^n5-point:([^:]+):([^:]+)$/.exec(id||'');
    return m?fiveQubitState(Number(m[1]),Number(m[2])):null;
  }
  function notify(text){$('notification').textContent=text;clearTimeout(notify.timer);notify.timer=setTimeout(()=>$('notification').textContent='',4500);}
  function showState(s,preview=false,{live=false}={}){if(!s)return;previewing=preview;document.querySelector('.majorana-section').hidden=false;$('no-state').hidden=true;$('state-formula').hidden=false;document.querySelector('.state-actions').hidden=false;$('numeric-details').hidden=false;current=stateData(s);const calc=sphere.setState(current.coefficients,{preserveOrientation:pointDragging});sector.setSelected(current);paintMath($('state-title'),AtlasMath.state(s.id,s.mathLabel||title(s)));$('selection-mode').textContent=preview?'Preview':'Selected';$('state-evidence').innerHTML=AtlasMath.html((s.proofOrigin||s.source||(s.status==='attained'?'Explicit normalized pure-state witness':'Explicit representative state')).replace(/\bS([1-6])\b/g,'$S_$1$'));showExpression(s,current.coefficients,{live});$('star-count').textContent=`${n} points · ${calc.groups.length} distinct`;const tuple=current.sectors;renderPurityValues(current);$('coefficients').innerHTML=`<table><thead><tr><th>k</th><th>c<sub>k</sub></th></tr></thead><tbody>${current.coefficients.map((c,k)=>`<tr><td>${k}</td><td>${esc(complex(c))}</td></tr>`).join('')}</tbody></table>`;const norm=current.coefficients.reduce((a,c)=>a+c[0]**2+c[1]**2,0),error=Math.max(...tuple.map((v,k)=>Math.abs(v-(s.sectors?.[k]??v))));$('verification').textContent=`Norm: ${fmt(norm,10)} · Sector check: ${error.toExponential(1)} · Constellation reconstruction fidelity: ${fmt(calc.fidelity,12)}${calc.reliable?'':' · Numerical constellation may be ill-conditioned'}.`;$('majorana-view').setAttribute('aria-label',`Majorana sphere for ${title(s)}: ${n} points, ${calc.groups.length} distinct positions. Drag to rotate.`);if(live)window.StatePicker?.updateLive?.(s);else{$('state-select').value=s.id||'';window.StatePicker?.update(s,n);}if(!preview&&!live)updateFamily(s);}
  function syncCompanionExplorer(state){
    const link=$('companion-explorer');
    link.hash=new URLSearchParams({n,state:state.id||''}).toString();
  }
  function pointRepresentatives(state){
    const dimensions=data.dimension;
    return named.filter(candidate=>candidate.sectors.slice(0,dimensions).length===dimensions &&
      candidate.sectors.slice(0,dimensions).every((v,k)=>Math.abs(v-state.sectors[k])<1e-8));
  }
  function updatePointRepresentatives(state){
    const panel=$('point-representatives'), options=$('point-representative-options');
    const matches=pointRepresentatives(state);
    panel.hidden=matches.length<2;
    if(panel.hidden)return;
    const ids=matches.map(s=>s.id).join(' ');
    if(options.dataset.ids!==ids){
      options.dataset.ids=ids;
      options.replaceChildren(...matches.map(candidate=>{
        const button=document.createElement('button');
        button.type='button';
        button.dataset.stateId=candidate.id;
        paintMath(button,AtlasMath.state(candidate.id,title(candidate)));
        button.addEventListener('click',()=>selectState(candidate));
        return button;
      }));
    }
    for(const button of options.children){
      button.setAttribute('aria-pressed',String(button.dataset.stateId===state.id));
    }
  }
  function selectState(s,updateURL=true,{live=false}={}){updatePointRepresentatives(s);syncCompanionExplorer(s);clearTimeout(previewTimer);previewing=false;if(live){selected=s;showState(s,false,{live:true});return;}if(s.interactivePoint){interactiveState=s;$('state-select').querySelectorAll('[data-interactive-point]').forEach(o=>o.remove());const op=document.createElement('option');op.value=s.id;op.dataset.interactivePoint='true';op.textContent=`Selected point: (${PurityCoordinates.fromSectors(n,s.sectors).map(v=>fmt(v,4)).join(', ')})`;$('state-select').append(op);}selected=s;showState(s);if(updateURL&&s.id){const h=`n=${n}&state=${encodeURIComponent(s.id)}`;history.replaceState(null,'','#'+h);}notify(`Selected ${title(s)}.`);}
  function updateFamily(s){$('family-control').hidden=true;}
  function familyText(s){const p=s.parameter||s.parameters||{};$('family-value').textContent=Object.entries(p).map(([k,v])=>`${k} = ${typeof v==='number'?fmt(v):v}`).join(' · ');}
  function sampleBoundary(d){const out=[];for(const curve of d.curves){const count=curve.points.length;for(let j=1;j<10;j++){const p=curve.points[Math.round(j*(count-1)/10)];if(!out.some(q=>Math.hypot(q.sectors[0]-p.sectors[0],q.sectors[1]-p.sectors[1])<.02)&&!d.knownStates.some(q=>Math.hypot(q.sectors[0]-p.sectors[0],q.sectors[1]-p.sectors[1])<.02))out.push({...p,curveId:curve.id,family:curve.label,label:p.label||curve.label,status:'attained',source:'Analytically attained boundary: '+curve.label});}}return out;}
  function getData(k){
    if(k!==2)return window.SECTOR_SMALL_DATA[String(k)];
    const knownStates=[
      {...intervalState(0,2),interactivePoint:false,id:'n2-ghz',exactSectors:['0','3'],label:'Bell state'},
      {...intervalState(2,2),interactivePoint:false,id:'n2-product',exactSectors:['2','1'],label:'Product state'}
    ];
    return {n:2,dimension:1,knownStates,curves:[]};
  }
  function setN(k,stateId){
    sector.cancelPointerGesture();
    pointDragging=false;
    n=Number(k);
    interactiveState=null;
    data=getData(n);
    named=(data.states||data.knownStates).filter(s=>s.id!=='n5-dicke2');
    if(n===5)data={...data,knownStates:named,boundarySamples:sampleBoundary(data)};
    document.querySelectorAll('[data-n]').forEach(b=>b.setAttribute('aria-pressed',Number(b.dataset.n)===n));
    $('geometry-title').textContent=names[n]+'-qubit geometry';
    $('set-status').textContent=data.dimension+' independent coordinate'+(data.dimension===1?'':'s')+' · Exact achievable set';
    const picker=$('state-select');
    picker.replaceChildren(new Option('Explore a point on the plot',''));
    for(const state of named)picker.append(new Option(title(state),state.id));
    if(data.boundarySamples){
      const group=document.createElement('optgroup');
      group.label='Sampled boundary states';
      for(const state of data.boundarySamples){
        const tuple=PurityCoordinates.fromSectors(n,state.sectors).map(v=>fmt(v,4)).join(', ');
        group.append(new Option(state.family+' ('+tuple+')',state.id));
      }
      picker.append(group);
    }
    $('sector-view').classList.toggle('line-view',data.dimension===1);
    $('attained-region-key').hidden=data.dimension===1;
    $('outer-label').textContent='Exact achievable set';
    $('plot-gesture').textContent='Drag the selected point · Drag elsewhere to move · Scroll to zoom';
    const note=n===2?'The full interval is $1/2\\le P_1\\le1$, with $P_2=1$.':n===3?'The full interval is $1/2\\le P_1\\le1$, with $P_2=P_1$.':
      n===4?'The full triangle is $1/2\\le P_1\\le1$ and $(4P_1-1)/3\\le P_2\\le P_1$.':
      'The achievable region is nonconvex.';
    $('geometry-note').innerHTML=AtlasMath.html(note);
    $('sector-canvas').setAttribute('aria-label',names[n]+'-qubit reduced-purity set. Click the set or drag the selected point to inspect an attaining state.');
    sector.setData(data);
    const initial=[...named,...(data.boundarySamples||[])].find(s=>s.id===stateId)||
      (n<=3?intervalStateFromId(stateId):n===4?triangleStateFromId(stateId):fiveQubitStateFromId(stateId))||
      named.find(s=>s.id==='n'+n+'-ghz')||named[0];
    selectState(initial,false);
    history.replaceState(null,'','#n='+n+'&state='+encodeURIComponent(initial.id));
    document.title='N = '+n+' · Reduced-purity geometry';
  }
  function preview(s){
    clearTimeout(previewTimer);
    if(pointDragging)return;
    if(!s){
      if(previewing)previewTimer=setTimeout(()=>{if(previewing&&selected)showState(selected);},100);
      return;
    }
    if(s.id&&s.id===selected?.id){if(previewing)showState(selected);return;}
    previewTimer=setTimeout(()=>showState(s,true),65);
  }
  function downloadJSON(record,filename){const blob=new Blob([JSON.stringify(record,null,2)],{type:'application/json'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=filename;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
  try{if(!window.THREE||!window.SectorMath||!window.SECTOR_SMALL_DATA||!window.PurityCoordinates)throw Error('A required local asset did not load. Keep the complete atlas folder together.');sector=new SectorView($('sector-canvas'),$('sector-labels'));sphere=new MajoranaView($('majorana-canvas'),$('majorana-labels'));function pinDisplayedState(){clearTimeout(previewTimer);if(previewing&&current){selected=current;updatePointRepresentatives(current);syncCompanionExplorer(current);previewing=false;$('selection-mode').textContent='Selected';if(current.id)history.replaceState(null,'',`#n=${n}&state=${encodeURIComponent(current.id)}`);}}document.querySelector('.majorana-section').addEventListener('pointerdown',pinDisplayedState,{capture:true});document.querySelector('.majorana-section').addEventListener('wheel',pinDisplayedState,{capture:true,passive:true});sector.hoverCallback=preview;sector.pointDragCallback=(coordinates,{phase})=>{if(phase==='start'){clearTimeout(previewTimer);pointDragging=true;previewing=false;if(current)selected=current;}else if(phase==='move'){const[x,y]=coordinates;const s=n<=3?intervalState(x):n===4?triangleState(x,y):n===5?fiveQubitState(x,y):null;if(s)selectState(s,false,{live:true});}else if(phase==='end'||phase==='cancel'){if(pointDragging&&selected)selectState(selected,true);pointDragging=false;}};sector.coordinateCallback=([x,y])=>{if(n>=2&&n<=5){const s=n<=3?intervalState(x):n===4?triangleState(x,y):fiveQubitState(x,y);if(s)selectState(s);}};sector.selectCallback=s=>{if(s&&!s.outerOnly){clearTimeout(previewTimer);selectState(s);}};document.querySelectorAll('[data-n]').forEach(b=>b.addEventListener('click',()=>setN(b.dataset.n)));document.querySelectorAll('[data-view]').forEach(b=>b.addEventListener('click',()=>{sector.setPreset(b.dataset.view);document.querySelectorAll('[data-view]').forEach(p=>p.setAttribute('aria-pressed',p===b));}));['show-outer','show-labels'].forEach(id=>$(id).addEventListener('change',()=>sector.applyOptions()));$('state-select').addEventListener('change',e=>{const s=[...named,...(data.boundarySamples||[]),...(interactiveState?[interactiveState]:[])].find(p=>p.id===e.target.value);if(s)selectState(s);});$('show-hull').addEventListener('change',()=>{sphere.hull.visible=$('show-hull').checked;sphere.render();});$('reset-sphere').addEventListener('click',()=>sphere.reset());document.querySelectorAll('[data-sphere-turn]').forEach(b=>b.addEventListener('click',()=>{const t=b.dataset.sphereTurn;if(t==='left')sphere.yaw-=.22;if(t==='right')sphere.yaw+=.22;if(t==='up')sphere.pitch=Math.min(1.3,sphere.pitch+.17);if(t==='down')sphere.pitch=Math.max(-1.3,sphere.pitch-.17);sphere.updateCamera();sphere.render();}));$('reset-sector-view').addEventListener('click',()=>{sector.resetView();document.querySelectorAll('[data-view]').forEach(b=>b.setAttribute('aria-pressed',b.dataset.view==='main'));});$('export-plot').addEventListener('click',()=>AtlasSVG.download(sector,{filename:`purity-set-N${n}.svg`,title:`Pure symmetric N=${n} reduced-purity set`}));$('export-majorana-svg').addEventListener('click',()=>AtlasSVG.download(sphere,{filename:`majorana-N${n}-${current.id||'state'}.svg`,title:`Majorana constellation: ${title(current)}`}));$('align-symmetry').addEventListener('click',()=>{clearTimeout(previewTimer);sphere.alignHighest();});$('show-symmetry').addEventListener('change',()=>sphere.render());$('symmetry-order').addEventListener('change',()=>sphere.render());$('download-state').addEventListener('click',()=>{const stars=SectorMath.majoranaStars(current.coefficients);downloadJSON({...current,normalization:'purities[k-1] is tr(rho_k^2). The auxiliary sectors array stores the sector lengths.',basis:'Normalized Dicke basis, k=0..N excitations; Jz=N/2-k.',majoranaConvention:SectorMath.convention||'Constituent spinor stars: |0> north, |1> south.',majorana:stars},`N${n}-${current.id||'selected-state'}.json`);});$('copy-link').addEventListener('click',async()=>{const url=new URL(location.href);if(current.id)url.hash=`n=${n}&state=${encodeURIComponent(current.id)}`;try{await navigator.clipboard.writeText(url.href);notify('Link copied.');}catch{notify('Copy the address from your browser to share this selection.');}});$('family-slider').addEventListener('input',e=>{if(!activeCurve)return;const s=activeCurve.points[Number(e.target.value)];selected=s;syncCompanionExplorer(s);showState(s);familyText(s);});window.addEventListener('hashchange',()=>{const q=new URLSearchParams(location.hash.slice(1)),k=Number(q.get('n'));if([2,3,4,5].includes(k))setN(k,q.get('state'));});const q=new URLSearchParams(location.hash.slice(1));setN([2,3,4,5].includes(Number(q.get('n')))?Number(q.get('n')):3,q.get('state'));AtlasMath.renderStatic();document.fonts.ready.then(()=>{sector.resize();sphere.resize();fitSectorValues();});window.addEventListener('resize',fitSectorValues);window.SectorAtlas={getN:()=>n,getSelected:()=>current,data:getData};}
  catch(e){console.error(e);const p=document.createElement('p');p.className='error';p.setAttribute('role','alert');p.textContent='The interactive viewer could not start. '+e.message;$('explorer').prepend(p);}
})();
