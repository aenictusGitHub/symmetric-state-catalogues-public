/* Finite 3D axes and label continuity, independent of the modal DOM harness. */
const fs=require('fs'),path=require('path');
const root=path.resolve(__dirname,'../../..');
const defaults=['anticoherent-state-catalogue/majorana-viewer.js','symmetric-sector-catalogue/majorana-viewer.js','ac-subspace-catalogue/majorana-viewer.js','majorana-degeneracy-catalogue/app.js'];
let assertions=0;
const assert=(ok,message)=>{assertions++;if(!ok)throw Error(message);};
const radius=180,cx=300,cy=225,width=600,height=450;
function drawAxis(source,angle,directions=new Map()){
 const labels=[],lines=[],markers=[];let current=[];
 const context={globalAlpha:1,setLineDash(){},save(){},restore(){},beginPath(){current=[];},moveTo(x,y){current.push([x,y]);},lineTo(x,y){current.push([x,y]);},stroke(){if(current.length===2)lines.push(current);},arc(x,y,r){markers.push({x,y,r,alpha:this.globalAlpha});}};
 const rotate=p=>({x:p.x,y:p.y*Math.cos(angle)-p.z*Math.sin(angle),z:p.y*Math.sin(angle)+p.z*Math.cos(angle)});
 const project=p=>({x:cx+radius*p.x,y:cy-radius*p.y,z:p.z});
 const axisLine=new Function('rotate','project','context','radius','cx','cy','width','height','axisLabels','axisDirections',source+';return axisLine;')(rotate,project,context,radius,cx,cy,width,height,labels,directions);
 axisLine([0,0,1],'C₈','#526b9b');
 assert(lines.length>=2&&lines.length<=4,'axis halves split only at sphere visibility boundaries');
 const distance=Math.max(...lines.flat().map(p=>Math.hypot(p[0]-cx,p[1]-cy)));
 assert(Math.abs(distance-radius*1.08*Math.abs(Math.sin(angle)))<1e-10,'axis foreshortens with its 3D projection');
 assert([...lines.flat(2),labels[0].x,labels[0].y].every(Number.isFinite),'finite draw coordinates');
 return {label:labels[0],distance,markers};
}
const reports=[];
for(const file of [...defaults,...process.argv.slice(2)]){
 const absolute=path.resolve(root,file),text=fs.readFileSync(absolute,'utf8'),start=text.indexOf('      function axisLine(vector,label,color,unoriented=true){'),end=text.indexOf('      if(controls.axes?.checked)',start);
 assert(start>=0&&end>start,'axis renderer found: '+file);
 const source=text.slice(start,end);let maxJump=0,previous=null;
 // Passing through the camera direction used to force a full-diameter
 // segment and a ~190px discontinuity at a four-pixel threshold.
 for(let i=-100;i<=100;i++){
  const value=drawAxis(source,i*.001);
  if(previous)maxJump=Math.max(maxJump,Math.hypot(value.label.x-previous.x,value.label.y-previous.y));
  previous=value.label;
 }
 assert(maxJump<1,'label stays continuous through an end-on axis');
 const before=drawAxis(source,Math.PI/2-.001),after=drawAxis(source,Math.PI/2+.001);
 assert(Math.hypot(before.label.x-after.label.x,before.label.y-after.label.y)<1,'label does not swap ends at front/back crossing');
 const near=drawAxis(source,.005),middle=drawAxis(source,.02),far=drawAxis(source,.05);
 assert(near.distance<middle.distance&&middle.distance<far.distance,'continuous finite-axis foreshortening');
 assert(near.markers[0].alpha>middle.markers[0].alpha,'end-on marker fades smoothly');
 reports.push({file,maximumLabelStep:maxJump,frontBackLabelJump:Math.hypot(before.label.x-after.label.x,before.label.y-after.label.y),passed:true});
}
console.log(JSON.stringify({passed:true,assertions,catalogues:reports},null,2));
