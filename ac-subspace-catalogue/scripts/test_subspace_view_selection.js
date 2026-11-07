/* Selection regressions against the public viewer and its real renderer.
 * Geometry fixtures exercise navigation, not catalogue-state existence claims.
 */
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const sourcePath=path.join(__dirname,'../subspace-viewer.js');
const source=fs.readFileSync(sourcePath,'utf8');
const escape=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fixture=(id,n,active)=>({id,n,dimension:2,t:1,label:id,active});
const first=fixture('selection-fixture-a',6,[4,6]),second=fixture('selection-fixture-b',8,[3,5,8]);
function createWorld(store=new Map()){
 let document,sequence=0;const frames=new Map(),timers=[],details=[];
 class Element{
  constructor(tag='div'){
   this.tagName=tag.toUpperCase();this.handlers=new Map();this.children=[];this.style={};this.checked=false;this.value='';this.disabled=false;this.hidden=false;this.textContent='';this.attributes={};this.dataset={};this.isConnected=true;this._html='';this.classes=new Set();
   this.classList={add:c=>this.classes.add(c),remove:c=>this.classes.delete(c),contains:c=>this.classes.has(c),toggle:(c,on)=>on?this.classes.add(c):this.classes.delete(c)};
  }
  set innerHTML(value){
   this._html=String(value);
   if(this.tagName==='DIALOG'){
    this.controls=new Map();
    for(const match of this._html.matchAll(/<([a-z]+)[^>]*\bid="([^"]+)"[^>]*>/g)){
     const element=new Element(match[1]);element.id=match[2];element.parentDialog=this;element.hidden=/\bhidden(?:[\s>])/.test(match[0]);element.checked=/\bchecked(?:[\s>])/.test(match[0]);
     for(const attr of match[0].matchAll(/([\w-]+)="([^"]*)"/g))element.setAttribute(attr[1],attr[2]);
     if(match[1]==='select'){
      const end=this._html.indexOf('</select>',match.index);element.innerHTML=this._html.slice(match.index+match[0].length,end);
     }
     this.controls.set(element.id,element);
    }
    const canvas=this.controls.get('sv-canvas');
    if(canvas){canvas.parentElement=new Element();canvas.parentElement.querySelector=()=>this.help||(this.help=new Element('span'));}
   }else if(this.tagName==='SELECT'){
    this.options=new Map();for(const match of this._html.matchAll(/<option\s+value="([^"]+)"[^>]*>([^<]*)/g)){const option=new Element('option');option.value=match[1];option.textContent=match[2];this.options.set(option.value,option);}this.value=this.options.keys().next().value||'';
   }else if(this.id==='sv-bars'){
    this.children=[];for(const match of this._html.matchAll(/<button[^>]*data-rank="([^"]+)"[^>]*>/g)){const button=new Element('button');button.setAttribute('data-rank',match[1]);button.disabled=/\bdisabled(?:[\s>])/.test(match[0]);this.children.push(button);}
   }
   if(this.tagName!=='DIALOG'&&this.parentDialog){
    for(const match of this._html.matchAll(/<([a-z]+)[^>]*\bid="([^"]+)"[^>]*>/g)){const element=new Element(match[1]);element.id=match[2];element.parentDialog=this.parentDialog;this.parentDialog.controls.set(element.id,element);}
   }

  }
  get innerHTML(){return this._html;}
  setAttribute(name,value){this.attributes[name]=String(value);if(name==='id')this.id=String(value);if(name.startsWith('data-'))this.dataset[name.slice(5)]=String(value);}
  getAttribute(name){return this.attributes[name]??null;}
  querySelector(selector){if(selector.startsWith('#'))return this.controls?.get(selector.slice(1))||null;const option=selector.match(/^option\[value="([^"]+)"\]$/);if(option)return this.options?.get(option[1])||null;return null;}
  querySelectorAll(selector){return selector==='button'?this.children.filter(e=>e.tagName==='BUTTON'):[];}
  addEventListener(name,fn){if(!this.handlers.has(name))this.handlers.set(name,new Set());this.handlers.get(name).add(fn);}
  removeEventListener(name,fn){this.handlers.get(name)?.delete(fn);}
  dispatch(name,extra={}){for(const fn of [...(this.handlers.get(name)||[])])fn({target:this,...extra});}
  click(){this.dispatch('click');}
  append(...children){this.children.push(...children);}
  replaceChildren(...children){this.children=[...children];}
  focus(){document.activeElement=this;}
  showModal(){this.open=true;this.querySelector('#sv-close')?.focus();}
  close(){this.open=false;this.dispatch('close');}
  getBoundingClientRect(){return {left:0,top:0,right:486,bottom:437,width:486,height:parseFloat(this.style.height)||437};}
  getContext(){return canvasContext;}
  closest(selector){if(selector==='dialog')return this.parentDialog||null;if(selector==='[data-rank]'&&this.getAttribute('data-rank')!=null)return this;return null;}
 }
 const canvasContext=new Proxy({measureText:s=>({width:String(s).length*7}),createRadialGradient:()=>({addColorStop(){}}),createImageData:(w,h)=>({data:new Uint8ClampedArray(w*h*4)})},{get:(o,k)=>k in o?o[k]:()=>{}});
 document={body:new Element('body'),activeElement:new Element('button'),createElement:tag=>new Element(tag)};
 const localStorage={getItem:key=>store.get(key)??null,setItem:(key,value)=>store.set(key,String(value))};
 const window={devicePixelRatio:1,localStorage,katex:{renderToString:tex=>`<span>${escape(tex)}</span>`},
  SubspaceSphereLayout:require('../subspace-sphere-layout.js'),HusimiRendering:require('../husimi-rendering.js'),SubspaceRenderGeometry:require('../subspace-render-geometry.js'),
  SubspaceSvgExport:{mathLabel:()=>''},SubspaceSphereDetail:{open:options=>details.push(options),close(){}},
  OperatorMath:{
   fromBasis:s=>({id:s.id,n:s.n,dimension:s.dimension,diagnostics:{},multipoles:Array.from({length:s.n+1},(_,rank)=>({rank,weight:s.active.includes(rank)?.1+rank/20:0,zero:!s.active.includes(rank),coefficients:Array.from({length:2*rank+1},()=>[0,0])}))}),
   husimi:model=>1/(model.n+1),symmetry:()=>({axes:[],tolerance:1e-6}),
   constellation:(model,rank)=>({rank,weight:model.multipoles[rank].weight,available:true,selectedPoints:[[1,0,0]],oppositePoints:[[-1,0,0]],rootResidual:0,reconstructionResidual:0,signReliable:true})
  },
  SpinWigner:{prepare:()=>({grid:()=>({rows:1,cols:1,qMin:-.1,qMax:.2,vertices:[{n:[1,0,0],q:-.1},{n:[0,1,0],q:0},{n:[0,0,1],q:.2},{n:[-1,0,0],q:0}],faces:[[0,1,2,3]]})})}
 };
 const scope={window,document,localStorage,console,Map,WeakMap,Date,Blob,URL,
  requestAnimationFrame:fn=>{frames.set(++sequence,fn);return sequence;},cancelAnimationFrame:id=>frames.delete(id),setTimeout:fn=>{timers.push(fn);return ++sequence;},ResizeObserver:class{observe(){}disconnect(){}}
 };
 vm.runInNewContext(source,scope,{filename:sourcePath});
 const dialog=()=>document.body.children.find(e=>e.id==='subspace-viewer');
 const ui=id=>dialog()?.querySelector('#sv-'+id);
 async function settle(){
  for(let i=0;i<8;i++){await Promise.resolve();const scheduled=[...frames.values()];frames.clear();scheduled.forEach(fn=>fn());timers.splice(0).forEach(fn=>fn());if(!frames.size&&!timers.length)await Promise.resolve();}
 }
 async function open(s){const promise=window.SubspaceViewer.open(s);await settle();await promise;await settle();}
 async function click(id){assert.ok(ui(id),`Missing ${id} control`);ui(id).click();await settle();}
 async function choose(rank){ui('rank').value=String(rank);ui('rank').dispatch('change');await settle();}
 const cells=()=>ui('math-labels').children.filter(e=>e.className==='sv-math-cell');
 const active=id=>ui(id).getAttribute('aria-pressed')==='true';
 const handlerCount=id=>[...ui(id).handlers.values()].reduce((sum,set)=>sum+set.size,0);
 return {window,document,store,dialog,ui,open,click,choose,settle,cells,active,handlerCount,details};
}
(async()=>{
 const store=new Map(),world=createWorld(store);await world.open(first);
 assert.equal(world.active('husimi'),true,'Fresh viewer should retain its Husimi default');
 assert.equal(world.ui('spheres').hidden,false,'Gallery entry must remain visible in the scalar view');
 await world.click('spheres');
 assert.equal(world.ui('rank').value,'all');assert.equal(world.ui('separate').checked,true);assert.equal(world.active('spheres'),true);assert.equal(world.active('operator'),false);
 assert.equal(world.cells().length,2,'One gallery cell is required for every nonzero multipole');
 assert.equal(store.get('ac-subspace-representation'),'spheres');
 const rankButton=world.cells()[0].children.find(e=>e.tagName==='BUTTON');assert.ok(rankButton,'Each gallery rank must remain clickable');rankButton.click();assert.equal(world.details.at(-1).rank,4);
 const listenersBefore=world.handlerCount('spheres');assert.equal(listenersBefore,1,'Gallery button must have exactly one action');
 world.window.SubspaceViewer.close();assert.equal(world.handlerCount('spheres'),0,'Closing must remove the prior gallery handler');
 await world.open(second);assert.equal(world.active('spheres'),true);assert.equal(world.ui('separate').checked,true);assert.equal(world.cells().length,3,'Changing subspace must rebuild the full gallery');assert.equal(world.handlerCount('spheres'),1,'Reopening must not accumulate gallery handlers');
 await world.click('operator');assert.equal(world.active('operator'),true);assert.equal(world.active('spheres'),false);assert.equal(world.ui('separate').checked,false,'Nested/individual view must remain available');assert.equal(world.cells().length,0);
 await world.choose(5);assert.equal(world.ui('rank').value,'5');assert.equal(world.ui('separate').checked,false);assert.equal(world.active('operator'),true);
 await world.click('husimi');assert.equal(world.active('husimi'),true);assert.equal(world.ui('spheres').hidden,false);
 await world.click('wigner');assert.equal(world.active('wigner'),true);assert.equal(world.ui('spheres').hidden,false);
 await world.click('spheres');assert.equal(world.cells().length,3);assert.equal(world.active('wigner'),false);
 // A fresh VM is a full page reload; only storage is carried between pages.
 const reloaded=createWorld(store);await reloaded.open(first);assert.equal(reloaded.active('spheres'),true,'A stored gallery preference must survive a page reload');assert.equal(reloaded.cells().length,2);
 await reloaded.click('husimi');reloaded.window.SubspaceViewer.close();await reloaded.open(second);assert.equal(reloaded.active('husimi'),true,'Explicitly switching away from the gallery must also persist');assert.equal(reloaded.cells().length,0);assert.equal(reloaded.ui('spheres').hidden,false);
 console.log('PASS: permanent gallery access, every nonzero multipole, clickable rank enlargement, distinct active states, subspace/page-reload persistence, individual/nested access and listener cleanup.');
})().catch(error=>{console.error(error);process.exitCode=1;});
