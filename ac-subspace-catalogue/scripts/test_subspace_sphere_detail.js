/* Lifecycle and camera-isolation checks against the real operator renderer.
 * A minimal DOM supplies native-dialog event ordering; no browser is controlled.
 * Stellar points here are geometry fixtures, not catalogue-state witnesses.
 */
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const escape=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
let document,clock=0;
const queuedClose=[],downloads=[],compositions=[];
class Element{
 constructor(tag='div'){
  this.tagName=tag.toUpperCase();this.handlers=new Map();this.children=[];this.style={};this.classes=new Set();this.checked=false;this.value='';this.disabled=false;this.hidden=false;this.textContent='';this._html='';this.attributes={};this.isConnected=true;
  this.classList={add:x=>this.classes.add(x),remove:x=>this.classes.delete(x),contains:x=>this.classes.has(x),toggle:(x,on)=>on?this.classes.add(x):this.classes.delete(x)};
 }
 set innerHTML(value){this._html=String(value);if(this.tagName==='DIALOG'){
   this.controls=new Map();for(const match of this._html.matchAll(/<([a-z]+)[^>]*\bid="([^"]+)"[^>]*>/g)){const e=new Element(match[1]);e.id=match[2];e.parentDialog=this;this.controls.set(e.id,e);}
   const canvas=this.controls.get('ssd-canvas');if(canvas){canvas.parentElement=new Element();canvas.parentElement.querySelector=()=>new Element('span');}
  }}
 get innerHTML(){return this._html;}
 setAttribute(k,v){this.attributes[k]=String(v);if(k==='id')this.id=v;}
 getAttribute(k){return this.attributes[k]??null;}
 querySelector(selector){return this.controls?.get(selector.replace(/^#/,''))||null;}
 addEventListener(name,fn){if(!this.handlers.has(name))this.handlers.set(name,new Set());this.handlers.get(name).add(fn);}
 removeEventListener(name,fn){this.handlers.get(name)?.delete(fn);}
 dispatch(name,event={}){for(const fn of [...(this.handlers.get(name)||[])])fn({target:this,...event});}
 append(...items){this.children.push(...items);}
 replaceChildren(...items){this.children=[...items];}
 getBoundingClientRect(){return {x:0,y:0,left:0,top:0,right:600,bottom:500,width:600,height:parseFloat(this.style.height)||500};}
 getContext(){return canvasContext;}
 closest(selector){return selector==='dialog'?this.parentDialog||null:null;}
 focus(){document.activeElement=this;this.focusCalls=(this.focusCalls||0)+1;}
 click(){if(this.tagName==='A')downloads.push({href:this.href,filename:this.download});this.dispatch('click');}
 showModal(){this.open=true;document.activeElement=this.querySelector('#ssd-close');}
 close(){this.open=false;queuedClose.push(()=>this.dispatch('close'));}
}
const canvasContext=new Proxy({measureText:s=>({width:String(s).length*7}),createRadialGradient:()=>({addColorStop(){}})},{get:(o,k)=>k in o?o[k]:()=>{}});
document={body:new Element('body'),activeElement:null,createElement:tag=>new Element(tag)};
const mainDialog=new Element('dialog');mainDialog.open=true;
const opener=new Element('button');opener.parentDialog=mainDialog;document.activeElement=opener;
const window={devicePixelRatio:1,katex:{renderToString:tex=>`<span>${escape(tex)}</span>`},
 SubspaceSphereLayout:require('../subspace-sphere-layout.js'),HusimiRendering:require('../husimi-rendering.js'),SubspaceRenderGeometry:require('../subspace-render-geometry.js'),
 SubspaceSvgExport:{mathLabel:(tex,x,y)=>`<text x="${x}" y="${y}">${escape(tex)}</text>`,compose:(view,options)=>{compositions.push(options);return view;}}
};
const scope={window,document,console,Map,Date,Blob,URL:{createObjectURL:()=>`blob:test-${++clock}`,revokeObjectURL(){}},requestAnimationFrame:()=>++clock,cancelAnimationFrame(){},ResizeObserver:class{observe(){}disconnect(){}},setTimeout(){}};
vm.createContext(scope);
const viewerPath=path.join(__dirname,'../subspace-viewer.js');
vm.runInContext(fs.readFileSync(viewerPath,'utf8'),scope,{filename:viewerPath});
assert.equal(typeof window.SubspaceSphereRenderer,'function','The reusable renderer must be exported and accept its own controls.');
const realFactory=window.SubspaceSphereRenderer;
const controls=new Map(),mainUi=id=>{if(!controls.has(id))controls.set(id,new Element(id==='canvas'?'canvas':'div'));return controls.get(id);};
mainUi('fit').checked=true;mainUi('shape').value='sphere';mainUi('order').value='highest';mainUi('canvas').parentElement=new Element();mainUi('canvas').parentElement.querySelector=()=>new Element('span');
const records=[{rank:4,weight:.3},{rank:8,weight:Math.sqrt(1365/13376)}].map(r=>({...r,available:true,selectedPoints:[[1,0,0],[0,1,0],[0,0,1]],oppositePoints:[[-1,0,0],[0,-1,0],[0,0,-1]]}));
const model={n:13,dimension:2,multipoles:records.map(r=>({rank:r.rank,weight:r.weight,zero:false}))};
const mesh={qMin:0,qMax:1,vertices:[{n:[1,0,0],q:0},{n:[0,1,0],q:.2},{n:[0,0,1],q:.5},{n:[-1,0,0],q:1}],faces:[[0,1,2,3]],contours:{segments:[]}};
const symmetry={axes:[{order:4,axis:[0,0,1]}]};
const main=realFactory(model,mesh,symmetry,mainUi);main.set('operator',4,records[0]);main.align();
const positions=renderer=>[...renderer.exportSvg().matchAll(/<circle cx="([^"]+)" cy="([^"]+)" r="5\.2"/g)].map(m=>m[1]+','+m[2]);
const mainBefore=positions(main);assert.ok(mainBefore.length>0,'Fixture must contain displayed stars.');
const instances=[];
window.SubspaceSphereRenderer=(...args)=>{const renderer=realFactory(...args),entry={renderer,ui:args[3],disposed:0};instances.push(entry);return {...renderer,dispose(){entry.disposed++;renderer.dispose();}};};
const detailPath=path.join(__dirname,'../subspace-sphere-detail.js');vm.runInContext(fs.readFileSync(detailPath,'utf8'),scope,{filename:detailPath});
let closed=0;
const open=extra=>window.SubspaceSphereDetail.open({model,mesh,symmetry,record:records[1],rank:8,weightTex:'\\sqrt{\\frac{1365}{13376}}',exactWeight:true,weightNumeric:records[1].weight,subspaceId:'n13-fixture',onClose:()=>closed++,...extra});
open();
const dialog=document.body.children.find(e=>e.id==='subspace-sphere-detail');assert.ok(dialog?.open);
assert.equal(instances.length,1);assert.notEqual(instances[0].ui('canvas'),mainUi('canvas'),'Enlarged view must own its canvas and controls.');
assert.ok(dialog.querySelector('#ssd-weight').innerHTML.includes('\\sqrt{\\frac{1365}{13376}}'),'Exact analytical expression must survive unchanged.');
assert.ok(dialog.querySelector('#ssd-title').innerHTML.includes('\\ell=8'));
assert.equal(dialog.querySelector('#ssd-fit').checked,true);
const before=positions(instances[0].renderer);dialog.querySelector('#ssd-right').dispatch('click');
assert.notDeepEqual(positions(instances[0].renderer),before,'The enlarged view should rotate.');assert.deepEqual(positions(main),mainBefore,'Rotating the enlarged view must preserve the main camera.');
dialog.querySelector('#ssd-svg').dispatch('click');assert.equal(downloads.at(-1).filename,'n13-fixture-multipole-8.svg');
assert.ok(compositions.at(-1).formulaRows[0].tex.includes('\\sqrt{\\frac{1365}{13376}}'),'SVG legend needs the same exact weight.');
assert.ok(compositions.at(-1).symmetryRows.some(r=>r.tex.includes('C_{4}')));
window.SubspaceSphereDetail.close();assert.equal(instances[0].disposed,1);assert.equal(closed,1);assert.equal(document.activeElement,opener);
// A native close event is queued. Opening another rank before it fires must
// not dispose the new renderer or restore focus out of its active dialog.
const smallWeight=1.2e-7,newModel={...model,multipoles:[{rank:4,weight:smallWeight,zero:false},model.multipoles[1]]};
open({model:newModel,rank:4,record:{...records[0],weight:smallWeight},weightTex:'1.2e-7',exactWeight:false,weightNumeric:smallWeight});
for(const callback of queuedClose.splice(0))callback();
assert.equal(dialog.open,true);assert.equal(instances[1].disposed,0);assert.equal(closed,1);
assert.ok(dialog.querySelector('#ssd-weight').innerHTML.includes('1.2\\times10^{-7}'),'Numerical scientific notation must remain mathematical.');
dialog.querySelector('#ssd-svg').dispatch('click');assert.ok(compositions.at(-1).formulaRows[0].tex.includes('1.2\\times10^{-7}'));
assert.equal(compositions.at(-1).formulaRows[0].note,'Numerical value');
assert.deepEqual(positions(main),mainBefore);
// Native Escape closes only this top dialog and dispatches its close event.
let stopped=false;dialog.dispatch('cancel',{stopPropagation(){stopped=true;}});assert.equal(stopped,true);
dialog.open=false;dialog.dispatch('close');assert.equal(instances[1].disposed,1);assert.equal(closed,2);assert.equal(mainDialog.open,true);
assert.deepEqual(positions(main),mainBefore);
// Selection changes can explicitly close the detail after its parent has
// closed. Do not restore focus into that closed parent or dispose twice.
mainDialog.open=true;document.activeElement=opener;open();mainDialog.open=false;
const focusCalls=opener.focusCalls||0;window.SubspaceSphereDetail.close();window.SubspaceSphereDetail.close();
for(const callback of queuedClose.splice(0))callback();
assert.equal(instances[2].disposed,1);assert.equal(closed,3);assert.equal(opener.focusCalls||0,focusCalls);
main.dispose();
console.log('PASS: actual independently scoped cameras, close/Escape/reopen lifecycle, parent-close cleanup, exact/numerical weights and SVG legends.');
