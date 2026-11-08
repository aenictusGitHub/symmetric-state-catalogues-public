/* Actual Canvas axis routine and SVG recorder: no browser or state-data mocks. */
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const root=path.resolve(__dirname,'..'),catalogues=['symmetric-sector-catalogue','anticoherent-state-catalogue','ac-subspace-catalogue'];
const radius=180,cx=300,cy=225,width=600,height=450,near=(a,b,e=1e-9)=>assert.ok(Math.abs(a-b)<=e,`${a} != ${b}`);
function sourceFor(name){
 const source=fs.readFileSync(path.join(root,name,'majorana-viewer.js'),'utf8'),start=source.indexOf('      function axisLine('),end=source.indexOf('      if(controls.axes?.checked)',start);
 assert.ok(start>=0&&end>start);return source.slice(start,end);
}
function draw(source,angle,{coordinate=false,context:provided}={}){
 const labels=[],strokes=[],markers=[],stack=[];let current=[];
 const context=provided||{globalAlpha:1,lineDashOffset:0,dash:[],setLineDash(v){this.dash=v;},save(){stack.push({alpha:this.globalAlpha,dash:this.dash,offset:this.lineDashOffset});},restore(){const s=stack.pop();this.globalAlpha=s.alpha;this.dash=s.dash;this.lineDashOffset=s.offset;},beginPath(){current=[];},moveTo(x,y){current.push([x,y]);},lineTo(x,y){current.push([x,y]);},stroke(){if(current.length===2)strokes.push({points:current,alpha:this.globalAlpha,dash:this.dash.slice(),offset:this.lineDashOffset});},arc(x,y,r){markers.push({x,y,r,alpha:this.globalAlpha});}};
 const rotate=p=>({x:p.x,y:p.y*Math.cos(angle)-p.z*Math.sin(angle),z:p.y*Math.sin(angle)+p.z*Math.cos(angle)}),project=p=>({x:cx+radius*p.x,y:cy-radius*p.y,z:p.z});
 const fn=new Function('rotate','project','context','radius','cx','cy','width','height','axisLabels','axisDirections',source+';return axisLine;')(rotate,project,context,radius,cx,cy,width,height,labels,new Map());
 fn([0,0,1],coordinate?'z':'C₄','#a45635',!coordinate,coordinate?null:'4');assert.ok(labels.every(label=>label.order===(coordinate?null:'4')),'Math order travels with the unchanged axis geometry');return {labels,strokes,markers};
}
const report=[];
for(const name of catalogues){
 const source=sourceFor(name);
 for(const angle of [.7,1.5,Math.PI/2+.03,2.4]){
  const {strokes,labels}=draw(source,angle),length=radius*Math.abs(Math.sin(angle));
  assert.ok(strokes.some(s=>Math.abs(s.alpha-.8*.22)<1e-12),'Interior sections fade');
  assert.ok(strokes.some(s=>s.alpha===.8),'Exposed sections stay clear');
  for(const stroke of strokes){
   const [a,b]=stroke.points,lo=Math.hypot(a[0]-cx,a[1]-cy)/length,hi=Math.hypot(b[0]-cx,b[1]-cy)/length,rear=stroke.dash.length>0;
   const boundary=rear?1/Math.abs(Math.sin(angle)):1;
   assert.ok(hi<=1.08+1e-10&&hi>lo);
   near(stroke.alpha,.8*((lo+hi)/2<boundary?.22:1));
   near(stroke.offset,-length*lo);
   if(lo>1e-10)near(lo,boundary);
   if(hi<1.08-1e-10)near(hi,boundary);
  }
  assert.ok(labels.every(p=>Number.isFinite(p.x)&&Number.isFinite(p.y)));
  const coordinates=draw(source,angle,{coordinate:true}).strokes;
  assert.equal(coordinates.length,2);assert.ok(coordinates.every(s=>s.alpha===.8&&s.offset===0),'Coordinate guides unchanged');
 }
 // Both exactly end-on rays and tiny nearby projections remain finite.
 for(const angle of [0,1e-10,Math.PI]){
  const result=draw(source,angle);assert.ok(result.strokes.flatMap(s=>s.points.flat()).every(Number.isFinite));
  assert.ok(result.markers.length);assert.ok(result.markers.every(p=>p.alpha<=.8*.22+1e-12));
 }
 let previous=null,maxLabelStep=0;
 for(let i=-100;i<=100;i++){
  const label=draw(source,i*.001).labels[0];if(previous)maxLabelStep=Math.max(maxLabelStep,Math.hypot(label.x-previous.x,label.y-previous.y));previous=label;
 }
 assert.ok(maxLabelStep<1,'Label is continuous through end-on direction');
 const before=draw(source,Math.PI/2-.001).labels[0],after=draw(source,Math.PI/2+.001).labels[0];assert.ok(Math.hypot(before.x-after.x,before.y-after.y)<1,'Labels do not swap at horizon');
 const window={};vm.runInNewContext(fs.readFileSync(path.join(root,name,'majorana-svg.js'),'utf8'),{window});
 const svg=new window.MajoranaSVG.Context();draw(source,1.5,{context:svg});const xml=svg.toSVG(width,height,'Axis','Unit sphere');
 assert.ok(xml.includes('opacity="0.176"')&&xml.includes('opacity="0.8"'),'SVG preserves faded and exposed sections');
 assert.ok(/stroke-dasharray="4 4" stroke-dashoffset="-/.test(xml),'SVG keeps rear dash phase');
 assert.ok(!/NaN|Infinity/.test(xml));
 report.push({catalogue:name,maximumLabelStep:maxLabelStep,SVGOpacityAndDashPhase:'PASS'});
}
console.log(JSON.stringify({status:'PASS',checks:['actual ray/sphere intersections','rear silhouette','coordinate guides unchanged','end-on marker fading','label continuity','Canvas/SVG opacity','dash phase'],catalogues:report}));
