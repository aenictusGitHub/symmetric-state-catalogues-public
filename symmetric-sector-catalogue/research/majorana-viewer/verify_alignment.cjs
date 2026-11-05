// Exercise real viewer controls and geometry without a browser dependency.
const fs=require('fs'),vm=require('vm'),path=require('path'),assert=require('assert/strict');
const root=process.argv[2]?path.resolve(process.argv[2]):path.resolve(__dirname,'../..');
let pending=null,labels=[],hullFills=0,sphereRadii=[];
const win={addEventListener(){},removeEventListener(){},devicePixelRatio:1};
const ctx={window:win,requestAnimationFrame:f=>(pending=f,1),cancelAnimationFrame(){pending=null;}};
vm.runInNewContext(fs.readFileSync(path.join(root,'majorana-hull.js'),'utf8'),ctx);
let source;
if(fs.existsSync(path.join(root,'majorana-viewer.js')))source=fs.readFileSync(path.join(root,'majorana-viewer.js'),'utf8');
else {
 const app=fs.readFileSync(path.join(root,'app.js'),'utf8');
 source='(() => { const symmetryColor=()=>"#34786d";'+app.slice(app.indexOf('  // Axes are unoriented.'),app.indexOf('  function renderProofs()'))+'  window.MajoranaViewer={};})();';
 assert.ok(app.includes("align:$('viewer-align'),hull:$('viewer-hull')"),'degeneracy controls connected');
}
vm.runInNewContext(source.replace('  window.MajoranaViewer=','  window.test={principalAxis,axisToZ,attachSphere};\n  window.MajoranaViewer='),ctx);
const {principalAxis,axisToZ,attachSphere}=win.test;
const apply=(q,v)=>{
 const [w,x,y,z]=q,[a,b,c]=v;
 return [(1-2*(y*y+z*z))*a+2*(x*y-w*z)*b+2*(x*z+w*y)*c,2*(x*y+w*z)*a+(1-2*(x*x+z*z))*b+2*(y*z-w*x)*c,2*(x*z-w*y)*a+2*(y*z+w*x)*b+(1-2*(x*x+y*y))*c];
};
const near=(a,b)=>assert.ok(Math.abs(a-b)<1e-10,`${a} != ${b}`);
for(const v of [[0,0,1],[0,0,-1],[1,0,0],[0,1,0],[2,-3,4],[-7,2,-5]]){
 const q=axisToZ(v),n=Math.hypot(...v),out=apply(q,v.map(x=>x/n));
 near(Math.hypot(...q),1);out.forEach((x,i)=>near(x,i===2?1:0));
 near(Math.hypot(...apply(q,[.2,.3,.4])),Math.hypot(.2,.3,.4));
}
const axes=[{order:2,direction:[0,0,1]},{order:3,direction:[1,2,-3].map(x=>x/Math.sqrt(14))},{order:3,direction:[1,0,0]}];
assert.equal(principalAxis({axes}).order,3);
assert.ok(principalAxis({axes}).direction[2]>0);
assert.equal(principalAxis({axes:[...axes,{order:null,direction:[1,0,0]}]}).order,null);
assert.equal(principalAxis({axes:[]}),null);
assert.equal(principalAxis({axes:[{order:3,direction:[NaN,0,1]}]}),null);
function el(){return {handlers:{},value:'100',checked:false,addEventListener(k,f){this.handlers[k]=f;},removeEventListener(k){delete this.handlers[k];}};}
const drawing=new Proxy({arc(x,y,r){if(x===300&&y===200&&r>100)sphereRadii.push(r);},fill(){if(String(this.fillStyle).startsWith('rgba(58,126,119,'))hullFills++;},fillText:(text,x,y)=>labels.push({text,x,y}),createRadialGradient:()=>({addColorStop(){}})},{get:(o,k)=>o[k]||(()=>{})});
const canvas={...el(),getContext:()=>drawing,getBoundingClientRect:()=>({width:600,height:400}),classList:{remove(){}},hasPointerCapture:()=>false};
const controls={hull:el(),reset:el(),align:el(),axes:el(),symmetry:{axes},symmetryToggle:el(),symmetryOrder:el()};
const points=[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(([x,y,z])=>({x:x/Math.sqrt(3),y:y/Math.sqrt(3),z:z/Math.sqrt(3),multiplicity:1})),before=JSON.stringify({points,axes});
const flush=()=>{labels=[];const f=pending;pending=null;f?.();};
controls.symmetryToggle.checked=true;controls.symmetryOrder.value='all';
const dispose=attachSphere(canvas,points,controls);flush();
assert.ok(sphereRadii.length);sphereRadii.forEach(r=>near(r,Math.min(600*.365,400*.405)*.95));assert.equal(controls.hull.disabled,false);assert.equal(controls.hull.checked,true);assert.equal(hullFills,4,'hull enabled on first frame');
assert.equal(controls.axes.checked,false);
assert.ok(!labels.some(p=>['x','y','z','1'].includes(p.text)),'no coordinate axes or single-point labels by default');
near(labels.find(p=>p.text==='C₃').x,300);assert.ok(labels.find(p=>p.text==='C₃').y<200,'highest-order axis aligned on first frame');
controls.hull.checked=false;controls.hull.handlers.change();flush();assert.equal(hullFills,4,'hull can be hidden');
controls.hull.checked=true;controls.hull.handlers.change();flush();assert.equal(hullFills,8,'hull can be shown again');
assert.equal(controls.align.disabled,false);controls.align.handlers.click();flush();
assert.equal(controls.axes.checked,false,'align preserves coordinate-axis visibility');assert.equal(controls.symmetryOrder.value,'3');controls.axes.checked=true;controls.axes.handlers.change();flush();
const z=labels.find(p=>p.text==='z'),c=labels.find(p=>p.text==='C₃');
near(z.x,300);near(c.x,300);assert.ok(z.y<c.y&&c.y<200,'z and principal axis upright, with distinct labels');
assert.equal(JSON.stringify({points,axes}),before,'alignment preserves saved data');
canvas.handlers.keydown({key:'ArrowLeft',preventDefault(){}});flush();assert.ok(Math.abs(labels.find(p=>p.text==='z').x-300)>1);controls.reset.handlers.click();flush();near(labels.find(p=>p.text==='z').x,300);sphereRadii.forEach(r=>near(r,Math.min(600*.365,400*.405)*.95));
dispose();assert.equal(controls.align.handlers.click,undefined);assert.equal(controls.hull.handlers.change,undefined);
controls.symmetry={axes:[]};const dispose2=attachSphere(canvas,points,controls);assert.equal(controls.align.disabled,true);dispose2();
console.log(root+' PASS: axis-to-z rotations, highest/continuous axis selection, actual alignment/reset controls, labels, cleanup, and unchanged source coordinates.');
