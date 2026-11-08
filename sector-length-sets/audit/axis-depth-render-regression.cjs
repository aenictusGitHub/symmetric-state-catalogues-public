/* Actual Three.js axis geometry and vector export checks, without WebGL. */
'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict');
const root=path.resolve(__dirname,'..'),T=require(path.join(root,'vendor/three.min.js'));
const controls=new Map(),element=()=>({checked:true,value:'all',style:{},options:[],innerHTML:'',append(o){this.options.push(o);},replaceChildren(){this.options=[];}});
const get=id=>{if(!controls.has(id))controls.set(id,element());return controls.get(id);};
const axes=[{direction:[.6,.4,Math.sqrt(.48)],order:4},{direction:[0,-1,0],order:2}];
const scope={THREE:T,MajoranaSymmetry:{analyze:()=>({axes})},AtlasMath:{html:s=>s},document:{getElementById:get,createElement:element,documentElement:{}},getComputedStyle:()=>({getPropertyValue:()=>'#faf9f5'}),XMLSerializer:class{},ATLAS_EXPORT_FONT_CSS:'.katex{}',console};scope.window=scope;vm.createContext(scope);
vm.runInContext(fs.readFileSync(path.join(root,'viewer.js'),'utf8'),scope);
vm.runInContext(fs.readFileSync(path.join(root,'svg-export.js'),'utf8'),scope);
const view=Object.create(scope.MajoranaView.prototype);Object.assign(view,{scene:new T.Scene(),content:new T.Group(),symmetryAxesGroup:new T.Group(),displayRotation:new T.Quaternion(),camera:new T.OrthographicCamera(-2,2,1.6,-1.6,.01,100),target:new T.Vector3(),width:640,height:512,labels:[],svg:{innerHTML:''},markerSvg:{innerHTML:''}});
view.scene.add(view.content);view.content.add(view.symmetryAxesGroup);view.camera.position.set(0,0,7);view.camera.lookAt(0,0,0);view.camera.updateMatrixWorld();
view.addLabel=()=>({el:element()});view.setupSymmetry([]);get('symmetry-order').value='all';
const near=(a,b,eps=2e-7)=>assert.ok(Math.abs(a-b)<eps,`${a} vs ${b}`);
const radii=object=>{const a=object.geometry.getAttribute('position');return [0,1].map(i=>Math.hypot(a.getX(i),a.getY(i),a.getZ(i)));};
function check(){view.updateSymmetry();let visible=0,hidden=0;
 for(const part of view.axisSegments){const depth=new T.Vector3(...part.direction).applyQuaternion(view.displayRotation).dot(view.camera.position.clone().normalize())*part.sign,front=depth>=0,boundary=Math.min(1.12,front?1:1/Math.sqrt(Math.max(1e-30,1-depth*depth)));
  for(const [object,base,occluded,dashed]of [[part.solid,.75,false,false],[part.dashed,.45,false,true],[part.occludedSolid,.75,true,false],[part.occludedDashed,.45,true,true]]){
   near(object.material.opacity,base*(occluded?.22:1));assert.equal(object.userData.axisOccluded,occluded);
   const on=(dashed?!front:front)&&(occluded?boundary>0:boundary<1.12-1e-10);assert.equal(object.visible,on);
   if(on){const [a,b]=radii(object);near(a,occluded?0:boundary);near(b,occluded?boundary:1.12);const distance=object.geometry.getAttribute('lineDistance');near(distance.getX(0),a);near(distance.getX(1),b);if(occluded)hidden++;else visible++;}
  }
 }
 assert.ok(hidden>0);return {hidden,visible};
}
const first=check();assert.ok(first.visible>0);
const svg=scope.AtlasSVG.serialize(view,{title:'Symmetry-axis depth regression'});assert.ok(svg.includes('opacity="0.165"')&&svg.includes('opacity="0.099"'));assert.ok(svg.includes('opacity="0.75"')&&svg.includes('opacity="0.45"'));assert.ok(/stroke-dashoffset="-[\d.]+"/.test(svg),'Rear exposed section preserves phase');assert.ok(!/NaN|Infinity|<image/.test(svg));
view.displayRotation.setFromAxisAngle(new T.Vector3(1,2,3).normalize(),.72);view.symmetryAxesGroup.quaternion.copy(view.displayRotation);view.camera.position.set(3,2,5);view.camera.lookAt(0,0,0);view.camera.updateMatrixWorld();check();
get('symmetry-order').value='4';view.updateSymmetry();assert.ok(view.axisSegments.filter(p=>p.order===2).every(p=>[p.solid,p.dashed,p.occludedSolid,p.occludedDashed].every(o=>!o.visible)));
get('show-symmetry').checked=false;view.updateSymmetry();assert.equal(view.symmetryAxesGroup.visible,false);
console.log('PASS: atlas actual Three.js axis depth, sphere/silhouette boundaries, rotation, filters, opacity, physical dash phase and vector SVG.');
