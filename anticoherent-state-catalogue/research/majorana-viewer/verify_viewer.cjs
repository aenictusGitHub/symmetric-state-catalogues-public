const fs=require('fs'),vm=require('vm'),path=require('path');
const catalogueRoot=path.resolve(__dirname,'../..'),base=path.dirname(catalogueRoot);
let assertions=0;const assert=(v,s)=>{assertions++;if(!v)throw Error(s);};
function session(name){
 const root=path.join(base,name),elements=new Map(),frames=new Map(),windowListeners=new Map();let frameId=0,doc,drawnFrames=0,observerCount=0;
 const fills=[],texts=[],arcs=[];
 const drawingContext=new Proxy({globalAlpha:1,createRadialGradient:()=>({addColorStop(){}}),measureText:t=>({width:String(t).length*8}),fill(){if(['#a45a34','#ddb392'].includes(this.fillStyle))fills.push({color:this.fillStyle,alpha:this.globalAlpha});},fillText:t=>texts.push(t),arc:(...a)=>arcs.push(a),clearRect(){drawnFrames++;}},{get:(o,k)=>k in o?o[k]:()=>{}});
 const classSet=()=>{const s=new Set();return {add:v=>s.add(v),remove:v=>s.delete(v),contains:v=>s.has(v)};};
 class El{
  constructor(tag){this.tagName=tag.toUpperCase();this._id='';this.listeners=new Map();this.value='';this.textContent='';this.classList=classSet();this.dataset={};this.attributes={};this.isConnected=true;this.disabled=false;this.open=false;this.capture=new Set();this._html='';this.owned=[];}
  set id(v){this._id=v;elements.set(v,this);}get id(){return this._id;}
  set innerHTML(v){for(const id of this.owned)elements.delete(id);this.owned=[];this._html=v;for(const m of v.matchAll(/<(\w+)[^>]*\bid="([^"]+)"[^>]*>/g)){const e=new El(m[1]);e.id=m[2];this.owned.push(e.id);const value=m[0].match(/\bvalue="([^"]*)"/);if(value)e.value=value[1];e.checked=/\bchecked\b/.test(m[0]);}}
  get innerHTML(){return this._html;}
  querySelector(s){if(s.startsWith('#'))return elements.get(s.slice(1))||null;return null;}
  querySelectorAll(){return [];}
  addEventListener(k,f){if(!this.listeners.has(k))this.listeners.set(k,new Set());this.listeners.get(k).add(f);}
  removeEventListener(k,f){this.listeners.get(k)?.delete(f);}
  fire(k,ev={}){for(const f of [...(this.listeners.get(k)||[])])f({target:this,preventDefault(){},...ev});}
  focus(){if(!this.disabled&&this.isConnected)doc.activeElement=this;}
  setAttribute(k,v){this.attributes[k]=v;}getContext(){return drawingContext;}
  showModal(){this.open=true;}close(){if(!this.open)return;this.open=false;this.fire('close');}
  setPointerCapture(id){this.capture.add(id);}hasPointerCapture(id){return this.capture.has(id);}releasePointerCapture(id){this.capture.delete(id);this.fire('lostpointercapture',{pointerId:id});}
  getBoundingClientRect(){return {width:600,height:450,left:0,top:0,right:600,bottom:450};}
  remove(){this.isConnected=false;elements.delete(this.id);}
 }
 const body=new El('body');body.append=()=>{};
 doc={activeElement:null,body,createElement:tag=>new El(tag),getElementById:id=>elements.get(id)||null};
 const win={katex:require(path.join(root,'vendor/katex/katex.min.js')),devicePixelRatio:1,addEventListener(k,f){if(!windowListeners.has(k))windowListeners.set(k,new Set());windowListeners.get(k).add(f);},removeEventListener(k,f){windowListeners.get(k)?.delete(f);}};
 class ResizeObserver {constructor(f){this.f=f;this.active=false;}observe(){this.active=true;observerCount++;}disconnect(){if(this.active)observerCount--;this.active=false;}}
 const ctx={window:win,document:doc,ResizeObserver,console,requestAnimationFrame:f=>{frames.set(++frameId,f);return frameId;},cancelAnimationFrame:id=>frames.delete(id)};
 for(const file of ['majorana-symmetry.js','symmetry-groups.js','point-groups.js','majorana-hull.js','majorana-svg.js','majorana-viewer.js'])vm.runInNewContext(fs.readFileSync(path.join(root,file),'utf8'),ctx,{filename:file});
 return {api:win.MajoranaViewer,win,elements,frames,El,doc,fills,texts,arcs,windowListeners,get observers(){return observerCount;},get drawnFrames(){return drawnFrames;},flush(){for(const [id,f] of [...frames]){frames.delete(id);f();}},setPoints:points=>({n:points.reduce((a,p)=>a+p.multiplicity,0),points})};
}
const summaries=[];
for(const name of [path.basename(catalogueRoot)]){
 const root=path.join(base,name),data=JSON.parse(fs.readFileSync(path.join(root,'data/majorana.json'))),records=Array.isArray(data.records)?data.records:Object.values(data.records);
 const env=session(name),button=new env.El('button');button.focus();
 assert(!env.api.open({n:2,points:[]}),name+' rejects empty');
 assert(!env.api.open({n:2,points:[{x:0,y:0,z:1,multiplicity:1}]}),name+' rejects weight mismatch');
 assert(!env.api.open({n:1,points:[{x:0,y:0,z:0,multiplicity:1}]}),name+' rejects non-unit');
 let count=0,finiteGroups=0,continuousGroups=0;const groupClassifications=[];
 for(const state of records){
  button.focus();assert(env.api.open(state),'opens '+state.id);const el=id=>env.elements.get('mv-'+id),dialog=env.elements.get('majorana-viewer');
  assert(dialog.open,'dialog opens');assert(el('download-svg')&&!el('download-svg').disabled,'SVG download available');assert(env.doc.activeElement===el('sphere'),'canvas focus');assert(env.doc.body.classList.contains('mv-open'),'scroll lock');
  assert(el('title').textContent.includes(state.n+' qubits'),'n title');assert(el('summary').textContent.includes(state.points.length+' distinct'),'positions count');const originalEvidence=state.evidence||'Numerical coordinates of the selected representative.';
  const displayedEvidence=el('evidence').textContent+' '+el('coordinate-note').textContent;
  for(const sentence of originalEvidence.split('. '))assert(displayedEvidence.includes(sentence),'evidence remains available');
  assert(!el('evidence').textContent.includes('Coordinates are a numerical rendering of the saved state.'),'rendering details removed from caption');
  assert(!el('evidence').textContent.includes('coincident finite roots are not inferred by clustering'),'clustering details removed from caption');
  assert(el('coordinate-note').hidden===!el('coordinate-note').textContent,'empty coordinate note hidden on representative changes');
  assert(el('groups').innerHTML.includes('Point group')&&el('groups').innerHTML.includes('Rotational subgroup'),'both group roles');assert(el('groups').innerHTML.includes('This representative only'),'group scope');assert(el('method-copy').innerHTML.includes('checked geometrically'),'numerical symmetry');
  assert(!el('zoom')&&!el('zoom-value'),'zoom slider removed');assert(el('axes').checked===false,'default axes');assert(!el('labels'),'single-point label control removed');
  const before=env.drawnFrames;env.flush();assert(env.drawnFrames===before+1,'draw frame');
  const detection=env.win.MajoranaSymmetry.detect(state.points);assert(detection.status!=='invalid','valid symmetry '+state.id);assert(detection.totalMultiplicity===state.n,'sym weight');if(detection.status==='continuous')continuousGroups++;else finiteGroups++;
  const proper=env.win.MajoranaGroups.classify(detection),point=env.win.MajoranaPointGroups.classify(state.points,detection,proper),shared=env.api.symmetryFor(state.points);
  assert(point.status==='classified'&&proper.status==='classified','both groups classified '+state.id);
  assert(shared===env.api.symmetryFor(state.points),'shared symmetry cache reused '+state.id);
  assert(shared.rotational.symbol===proper.symbol&&shared.point.symbol===point.symbol,'shared groups agree with direct geometry classification '+state.id);
  const summary=env.api.groupSummary(state),modal=el('groups').innerHTML;
  for(const group of [point,proper]){const rendered=env.win.katex.renderToString(group.symbol,{throwOnError:false,trust:false,strict:'ignore'});assert(summary.includes(rendered),'summary contains expected group '+state.id+':'+group.symbol);assert(modal.includes(rendered),'modal contains same group '+state.id+':'+group.symbol);}
  assert(summary.includes('Numerical geometry')&&summary.includes('This representative'),'summary numerical evidence scope '+state.id);
  assert(env.api.groupLabel(state).includes(point.label),'plain option group matches full point group '+state.id);
  groupClassifications.push({id:state.id,pointGroup:point.label,rotationalSubgroup:proper.label,pointGroupOrder:point.order,rotationalOrder:proper.order,numerical:true});
  if(detection.axes.length>8)assert(el('symmetry-order').value!=='all','dense axes filtered');
  el('sphere').fire('keydown',{key:'ArrowLeft'});env.flush();el('right').fire('click');env.flush();
  el('sphere').fire('pointerdown',{button:0,pointerType:'mouse',pointerId:count,clientX:0,clientY:0});assert(el('sphere').hasPointerCapture(count),'capture acquired');
  el('sphere').fire('pointermove',{pointerId:count,clientX:20,clientY:15});env.flush();
  env.api.close();assert(!dialog.open,'modal close');assert(!env.doc.body.classList.contains('mv-open'),'unlock');assert(env.doc.activeElement===button,'restore focus');assert(!el('sphere').hasPointerCapture(count),'capture released');assert(!el('sphere').classList.contains('dragging'),'drag removed');assert(!env.frames.size,'no pending frames');assert(env.observers===0,'observer removed');assert(!env.windowListeners.get('resize')?.size,'window resize removed');assert([...el('sphere').listeners.values()].every(x=>x.size===0),'canvas listeners removed');for(const key of ['download-svg','reset','left','symmetry','symmetry-order','axes'])assert([...el(key).listeners.values()].every(x=>x.size===0),'control listeners removed '+key);
  count++;
 }
 // Replace an open view and exercise Home reset, showing no duplicate controls.
 button.focus();env.api.open(records[0]);env.api.open(records[records.length-1]);env.flush();let e=id=>env.elements.get('mv-'+id);assert(e('title').textContent.includes(records.at(-1).n+' qubits'),'open replacement');assert(e('reset').listeners.get('click').size===1,'one reset listener');e('sphere').fire('keydown',{key:'Home'});env.flush();assert(env.arcs.some(a=>Math.abs(a[2]-Math.min(600*.365,450*.405)*.95)<1e-10),'fixed 95% sphere radius');env.api.close();assert(env.doc.activeElement===button,'replacement retains original opener');
 assert(env.fills.some(x=>x.color==='#ddb392'&&x.alpha===.85),'rear fill');assert(env.fills.some(x=>x.color==='#a45a34'&&x.alpha===1),'front fill');
 summaries.push({catalogue:name,recordsOpened:count,drawnFrames:env.drawnFrames,finiteGroups,continuousGroups,rearPointFillVerified:true,allCleanupChecksPassed:true,groupClassifications});
}
const report={passed:true,assertions,catalogues:summaries,checks:['every representative opens with the matching count, caption and evidence','weighted symmetry inputs valid','both full and proper group labels displayed with numerical scope','drawing with stronger rear-star fill','rotation keys/buttons and drag','fixed 95% scale and Home reset','focus and body lock','repeated open replacement','all renderer listeners, resize observers, frames and pointer captures disposed','every compact/card group agrees with the modal and direct geometric classification','shared symmetry classification cache reused','plain group labels match point group for AC options','SVG download button enabled and its listener disposed'],limitation:'Minimal DOM/canvas harness. Native modal keyboard focus and visual layout must be checked by the real browser.'};
fs.writeFileSync(path.join(__dirname,'viewer-audit.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report,null,2));
