/* Real font-outline, renderer and standalone SVG regressions. No browser required. */
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const root=path.resolve(__dirname,'..'),folders=['symmetric-sector-catalogue','anticoherent-state-catalogue','ac-subspace-catalogue','majorana-degeneracy-catalogue'];
const near=(a,b,eps=1e-9)=>assert.ok(Math.abs(a-b)<=eps,`${a} != ${b}`);
class Outline {constructor(d){this.d=d;}}
function element(){const handlers=new Map();return {checked:false,value:'all',handlers,classList:{add(){},remove(){}},addEventListener:(name,fn)=>handlers.set(name,fn),removeEventListener:(name,fn)=>{if(handlers.get(name)===fn)handlers.delete(name);},dispatch(name,event={}){handlers.get(name)?.({preventDefault(){},...event});},getBoundingClientRect:()=>({width:500,height:440}),focus(){}};}
function canvasContext(){
 let state={globalAlpha:1,lineWidth:1,strokeStyle:'#000',fillStyle:'#000',lineDashOffset:0,dash:[]},stack=[],current=[];const records=[],glyphs=[],texts=[],textPlacements=[];
 const ctx={records,glyphs,texts,textPlacements,save(){stack.push({...state,dash:[...state.dash]});},restore(){state=stack.pop();},setLineDash(dash){state.dash=[...dash];},beginPath(){current=[];},moveTo(x,y){current.push(['M',x,y]);},lineTo(x,y){current.push(['L',x,y]);},arc(...args){current.push(['A',...args]);},stroke(){records.push({...state,dash:[...state.dash],path:current.map(x=>[...x])});},createRadialGradient(){return {addColorStop(){}};},fill(outline){if(outline instanceof Outline)glyphs.push({d:outline.d,color:state.fillStyle,alpha:state.globalAlpha});},fillText(value,x,y){texts.push(String(value));textPlacements.push({text:String(value),x,y,color:state.fillStyle,alpha:state.globalAlpha});},strokeText(value){texts.push(String(value));}};
 for(const k of ['globalAlpha','lineWidth','strokeStyle','fillStyle','lineDashOffset','font','textAlign','textBaseline','lineJoin','shadowColor','shadowBlur','shadowOffsetX','shadowOffsetY'])Object.defineProperty(ctx,k,{get:()=>state[k],set:v=>state[k]=v});
 for(const name of ['setTransform','clearRect','translate','scale','closePath'])ctx[name]=()=>{};
 return ctx;
}
function environment(){
 const elements=new Map(),window={devicePixelRatio:1,MAJORANA_CATALOGUE:{entries:[],states:[],numbers:[]},addEventListener(){},removeEventListener(){}};
 const document={getElementById:id=>{if(!elements.has(id))elements.set(id,element());return elements.get(id);},querySelectorAll:()=>[]};let sequence=0;const frames=new Map();
 const scope={window,document,Path2D:Outline,requestAnimationFrame:fn=>{frames.set(++sequence,fn);return sequence;},cancelAnimationFrame:id=>frames.delete(id),ResizeObserver:class{observe(){}disconnect(){}},setTimeout(){},console};window.Path2D=Outline;
 const flush=()=>{const pending=[...frames.values()];frames.clear();for(const fn of pending)fn();};vm.createContext(scope);return {scope,window,flush};
}
const results=[];
for(const folder of folders){
 const location=path.join(root,folder),{scope,window,flush}=environment();
 vm.runInContext(fs.readFileSync(path.join(location,'majorana-latex.js'),'utf8'),scope);
 const api=window.MajoranaLatex;assert.ok(api&&typeof api.draw==='function'&&typeof api.toSVG==='function'&&typeof api.layout==='function');
 assert.equal(api.metrics.unitsPerEm,1000);near(api.metrics.scriptScale,.7);near(api.metrics.subscriptShift,.15);
 assert.equal(api.glyphs.C.advance,715,'C uses the vendored KaTeX math-italic advance');for(const c of '0123456789')assert.equal(api.glyphs[c].advance,500,'Subscript digits use KaTeX Main metrics');assert.equal(api.glyphs['∞'].advance,1000);
 for(const [letter,advance]of Object.entries({x:572,y:490,z:465})){assert.equal(api.glyphs[letter].advance,advance,'Coordinates use the vendored KaTeX math-italic metrics');assert.equal(api.layout(letter).tex,letter);near(api.layout(letter).width,16*advance/1000);}
 const labels=['2','12','100','continuous','x','y','z'];
 for(const order of labels){
  const layout=api.layout(order,16),svg=api.toSVG(order,200,100,16,'#a45635'),ctx=canvasContext();api.draw(ctx,order,200,100,16,'#a45635');
  assert.ok(layout.width>0&&layout.height>0&&[layout.width,layout.height,layout.top,layout.bottom].every(Number.isFinite));
  assert.ok(svg.includes('data-tex='),'SVG carries the mathematical label');assert.ok(!/<image|foreignObject|<text|NaN|Infinity/.test(svg),'Glyph export is native finite vector geometry');
  if(/^[xyz]$/.test(order))assert.equal(ctx.glyphs.length,1,'Coordinate labels are single math-italic outlines');else assert.ok(ctx.glyphs.length>=2,'Canvas draws the italic C and subscript outlines');assert.equal((svg.match(/<path\b/g)||[]).length,ctx.glyphs.length,'Canvas and SVG draw identical glyph counts');
  for(const glyph of ctx.glyphs){assert.ok(glyph.d.length>20&&svg.includes(glyph.d),'Canvas and SVG use the same actual font outline');assert.equal(glyph.color,'#a45635');near(glyph.alpha,1);}
  assert.equal(ctx.texts.length,0,'Mathematical labels do not use Unicode text fallback');
 }
 const one=api.layout('2',16),two=api.layout('12',16),three=api.layout('100',16);near(two.width-one.width,16*.7*.5);near(three.width-two.width,16*.7*.5);
 assert.equal(api.layout('continuous').tex,api.layout('infinity').tex);assert.equal(api.layout('continuous').tex,api.layout('∞').tex);
 vm.runInContext(fs.readFileSync(path.join(location,'majorana-svg.js'),'utf8'),scope);
 for(const order of labels){const recorder=new window.MajoranaSVG.Context();recorder.translate(7,11);recorder.globalAlpha=.37;api.draw(recorder,order,200,100,16,'#a45635');const xml=recorder.toSVG(500,440,'LaTeX axis','Vector glyph test');assert.ok(xml.includes(`data-tex="${api.layout(order).tex}"`));assert.ok(xml.includes('transform="matrix(1 0 0 1 7 11)" opacity="0.37"'),'SVG preserves the canvas transformation and opacity');assert.ok(!/<image|foreignObject|<text/.test(xml));assert.equal((xml.match(/<path\b/g)||[]).length,canvasGlyphCount(order));}
 function canvasGlyphCount(order){const ctx=canvasContext();api.draw(ctx,order,200,100,16,'#a45635');return ctx.glyphs.length;}
 // Exercise each complete current renderer rather than a duplicate draw routine.
 let source=fs.readFileSync(path.join(location,folder==='majorana-degeneracy-catalogue'?'app.js':'majorana-viewer.js'),'utf8');
 if(folder==='majorana-degeneracy-catalogue'){
  const hook="  const symmetryColor=order=>({2:'#34786d',3:'#82609b',4:'#a45635',5:'#b38b23',continuous:'#356f94'}[order]||'#526b9b');";assert.equal(source.split(hook).length,2);source=source.replace(hook,hook+'\n window.AxisLatexTestAttach=attachSphere;return;');
 }else{const hook='  window.MajoranaViewer=Object.freeze(';assert.equal(source.split(hook).length,2);source=source.replace(hook,'  window.AxisLatexTestAttach=attachSphere;\n'+hook);}
 vm.runInContext(source,scope);
 const drawn=[],instrumented={...api,draw(...args){drawn.push({order:args[1],x:args[2],y:args[3],size:args[4],color:args[5]});return api.draw(...args);}};window.MajoranaLatex=instrumented;
 const ctx=canvasContext(),canvas=element();canvas.getContext=()=>ctx;const download=element(),toggle=element();toggle.checked=true;
 const controls={symmetry:{axes:[{direction:[0,0,1],order:2},{direction:[1,0,0],order:12},{direction:[0,1,0],order:100},{direction:[.6,0,.8],order:null}]},symmetryToggle:toggle,symmetryOrder:element(),axes:element(),download,reset:element(),turns:[]};
 let exported='';window.MajoranaSVG={...window.MajoranaSVG,download:svg=>exported=svg};
 const stop=window.AxisLatexTestAttach(canvas,[{x:0,y:0,z:1,multiplicity:1},{x:0,y:0,z:-1,multiplicity:1}],controls);flush();
 assert.equal(drawn.length,4,'Every shown symmetry axis uses the mathematical renderer');assert.ok(ctx.glyphs.length>=9);assert.ok(!ctx.texts.some(t=>/^C[∞₀₁₂₃₄₅₆₇₈₉]/.test(t)),'No legacy Unicode axis label is painted');
 const rays=ctx.records.filter(r=>r.lineWidth===1.8&&r.path.length===2);assert.ok(rays.length);assert.ok(rays.some(r=>Math.abs(r.globalAlpha-.8*.22)<1e-12));assert.ok(rays.every(r=>Math.abs(r.globalAlpha-.8)<1e-12||Math.abs(r.globalAlpha-.8*.22)<1e-12),'Glyph drawing preserves axis depth opacity');assert.ok(rays.some(r=>r.dash.length&&r.lineDashOffset<0));
 download.dispatch('click');assert.ok(exported.includes('data-tex='));assert.ok(!/<image|foreignObject|NaN/.test(exported));assert.equal((exported.match(/data-tex=/g)||[]).length,4,'Complete exports retain all mathematical axis labels');assert.ok(!/>C[∞₀₁₂₃₄₅₆₇₈₉]/.test(exported));
 controls.axes.checked=true;ctx.texts.length=0;ctx.records.length=0;drawn.length=0;controls.axes.dispatch('change');flush();
 assert.equal(drawn.length,7,'Coordinate and symmetry labels all use the mathematical renderer');assert.ok(!ctx.texts.some(t=>/^[xyz]$/.test(t)),'Coordinate labels never use plain canvas text when LaTeX is available');
 const coordinateCalls=drawn.filter(call=>/^[xyz]$/.test(call.order));assert.equal(coordinateCalls.length,3);
 for(const [letter,color]of Object.entries({x:'#a15435',y:'#41776e',z:'#76648a'})){const call=coordinateCalls.find(call=>call.order===letter);assert.ok(call);assert.equal(call.color,color);assert.equal(call.size,16);assert.ok(Number.isFinite(call.x)&&Number.isFinite(call.y));}
 const guides=ctx.records.filter(r=>r.lineWidth===1&&r.path.length===2&&['#a15435','#41776e','#76648a'].includes(r.strokeStyle));assert.equal(guides.length,6);assert.ok(guides.every(r=>r.globalAlpha===.8),'Coordinate guide transparency remains unchanged');
 // The supported no-helper fallback independently records the original anchors.
 // Compare those anchors with the new outlines without duplicating camera math.
 window.MajoranaLatex=undefined;ctx.textPlacements.length=0;controls.axes.dispatch('change');flush();
 for(const call of coordinateCalls){const old=ctx.textPlacements.find(t=>t.text===call.order);assert.ok(old);near(call.x,old.x);near(call.y,old.y);assert.equal(call.color,old.color);}
 window.MajoranaLatex=instrumented;download.dispatch('click');assert.equal((exported.match(/data-tex=/g)||[]).length,7);assert.ok(!/<image|foreignObject|NaN/.test(exported));
 for(const call of coordinateCalls){assert.ok(exported.includes(`data-tex="${call.order}"`));assert.ok(exported.includes(`fill="${call.color}"`));assert.ok(!new RegExp(`>${call.order}</text>`).test(exported),'Coordinate export uses native LaTeX paths');}
 stop();results.push({catalogue:folder,actualRenderer:'PASS',fontOutlines:'PASS',CanvasSVGParity:'PASS',coordinateAnchors:'PASS',axisDepthOpacity:'PASS'});
}
console.log(JSON.stringify({status:'PASS',labels:['C_2','C_12','C_100','C_infinity','x','y','z'],results}));
