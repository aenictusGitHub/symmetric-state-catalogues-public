/* Actual degeneracy renderer and SVG-recorder regression, no browser. */
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const root=path.resolve(__dirname,'../majorana-degeneracy-catalogue');
const near=(a,b,e=1e-8)=>assert.ok(Math.abs(a-b)<e,`${a} != ${b}`);
function element(){const handlers=new Map();return {checked:false,value:'all',handlers,classList:{add(){},remove(){}},addEventListener:(name,fn)=>handlers.set(name,fn),removeEventListener:(name,fn)=>{if(handlers.get(name)===fn)handlers.delete(name);},dispatch(name,event={}){handlers.get(name)?.({preventDefault(){},...event});},getBoundingClientRect:()=>({width:500,height:440}),focus(){}};}
function canvasContext(){
 let state={globalAlpha:1,lineWidth:1,strokeStyle:'#000',fillStyle:'#000',lineDashOffset:0,dash:[]},stack=[],current=[];const records=[];
 const ctx={records,save(){stack.push({...state,dash:[...state.dash]});},restore(){state=stack.pop();},setLineDash(dash){state.dash=[...dash];},beginPath(){current=[];},moveTo(x,y){current.push(['M',x,y]);},lineTo(x,y){current.push(['L',x,y]);},arc(...args){current.push(['A',...args]);},stroke(){records.push({...state,dash:[...state.dash],path:current.map(x=>[...x])});},createRadialGradient(){return {addColorStop(){}};}};
 for(const k of ['globalAlpha','lineWidth','strokeStyle','fillStyle','lineDashOffset','font','textAlign','textBaseline','lineJoin','shadowColor','shadowBlur','shadowOffsetX','shadowOffsetY'])Object.defineProperty(ctx,k,{get:()=>state[k],set:v=>state[k]=v});
 for(const name of ['setTransform','clearRect','translate','scale','fill','closePath','fillText','strokeText'])ctx[name]=()=>{};
 return ctx;
}
const window={MAJORANA_CATALOGUE:{entries:[],states:[],numbers:[]},addEventListener(){},removeEventListener(){}},elements=new Map();
const document={getElementById:id=>{if(!elements.has(id))elements.set(id,element());return elements.get(id);},querySelectorAll:()=>[]};
let nextFrame=0;const frames=new Map();
const environment={window,document,requestAnimationFrame:fn=>{frames.set(++nextFrame,fn);return nextFrame;},cancelAnimationFrame:id=>frames.delete(id),ResizeObserver:class{observe(){}disconnect(){}},setTimeout(){},Map,Set,console};
let source=fs.readFileSync(path.join(root,'app.js'),'utf8');const hook="  const symmetryColor=order=>({2:'#34786d',3:'#82609b',4:'#a45635',5:'#b38b23',continuous:'#356f94'}[order]||'#526b9b');";
assert.equal(source.split(hook).length,2);source=source.replace(hook,hook+'\n  window.TestAttachSphere=attachSphere;return;');vm.runInNewContext(source,environment);
vm.runInNewContext(fs.readFileSync(path.join(root,'majorana-svg.js'),'utf8'),environment);
let exported='';window.MajoranaSVG={...window.MajoranaSVG,download:svg=>exported=svg};
const flush=()=>{const pending=[...frames.values()];frames.clear();for(const fn of pending)fn();};
const ctx=canvasContext(),canvas=element();canvas.getContext=()=>ctx;const controls={symmetry:{axes:[{direction:[0,0,1],order:4}]},symmetryToggle:element(),symmetryOrder:element(),axes:element(),download:element(),reset:element(),turns:[]};controls.symmetryToggle.checked=true;
const points=[{x:0,y:0,z:1,multiplicity:1},{x:0,y:0,z:-1,multiplicity:1}];
const dispose=window.TestAttachSphere(canvas,points,controls);flush();
const axisLines=()=>ctx.records.filter(r=>r.strokeStyle==='#a45635'&&r.lineWidth===1.8&&r.path.length===2&&r.path.every(p=>p[0]==='M'||p[0]==='L'));
let lines=axisLines();assert.equal(lines.length,4,'Both front and rear rays have hidden/exposed sections at the initial tilt');
const radius=Math.min(500*.365,440*.405)*.95,cx=250,cy=220,length=radius*Math.cos(.28),distance=p=>Math.hypot(p[1]-cx,p[2]-cy);
for(const r of lines){assert.ok(r.path.flat().slice(1).filter(x=>typeof x==='number').every(Number.isFinite));if(r.globalAlpha<.5)near(r.globalAlpha,.8*.22);else near(r.globalAlpha,.8);}
const front=lines.filter(r=>!r.dash.length),rear=lines.filter(r=>r.dash.length);assert.equal(front.length,2);assert.equal(rear.length,2);
near(distance(front.find(r=>r.globalAlpha<.5).path[1]),length);near(distance(front.find(r=>r.globalAlpha>.5).path[0]),length);
near(distance(rear.find(r=>r.globalAlpha>.5).path[0]),radius);assert.ok(rear.find(r=>r.globalAlpha>.5).lineDashOffset<0,'Split rear stroke retains a nonzero continuous dash phase');
controls.download.dispatch('click');assert.ok(exported.includes('opacity="0.176"'));assert.ok(exported.includes('opacity="0.8"'));assert.ok(/stroke-dashoffset="-/.test(exported));assert.ok(!/NaN|Infinity/.test(exported));
const svgAxis=[...exported.matchAll(/<path[^>]*stroke="#a45635"[^>]*>/g)].map(x=>x[0]);assert.equal(svgAxis.length,lines.length);for(const a of [.176,.8])assert.equal(svgAxis.filter(s=>s.includes(`opacity="${a}"`)).length,2);
// Turning changes the depth test; the same source axis does not teleport.
ctx.records.length=0;for(let i=0;i<9;i++){canvas.dispatch('keydown',{key:'ArrowDown'});flush();}lines=axisLines();assert.ok(lines.some(r=>r.globalAlpha<.5));assert.ok(lines.some(r=>r.globalAlpha>.5));assert.ok(lines.every(r=>r.path.every(p=>p.slice(1).every(Number.isFinite))));
// Coordinate guides do not acquire the new attenuation, even inside the globe.
controls.symmetryToggle.checked=false;controls.symmetryToggle.dispatch('change');controls.axes.checked=true;ctx.records.length=0;controls.axes.dispatch('change');flush();assert.equal(axisLines().length,0);const coordinate=ctx.records.filter(r=>r.lineWidth===1&&['#a15435','#41776e','#76648a'].includes(r.strokeStyle));assert.equal(coordinate.length,6);assert.ok(coordinate.every(r=>r.globalAlpha===.8));
// The continuous transverse C2 family uses a single unit-sphere ring.
const familyCtx=canvasContext(),familyCanvas=element();familyCanvas.getContext=()=>familyCtx;const familyToggle=element();familyToggle.checked=true;
const stopFamily=window.TestAttachSphere(familyCanvas,points,{symmetry:{axes:[{direction:[0,0,1],order:null}],continuousFamilies:[{}]},symmetryToggle:familyToggle,symmetryOrder:element()});flush();
const ring=familyCtx.records.filter(r=>r.strokeStyle==='#34786d'&&r.lineWidth===2.2);assert.equal(ring.length,180);assert.ok(ring.some(r=>r.globalAlpha===.22));assert.ok(ring.some(r=>r.globalAlpha===1));
stopFamily();dispose();
console.log('PASS: actual degeneracy inline/dialog sphere axes, exact front/rear emergence, opacity, rotation, dash phase, coordinate preservation, continuous family and SVG parity.');
