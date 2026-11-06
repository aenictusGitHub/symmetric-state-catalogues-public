'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const root=path.resolve(__dirname,'..'),T=require(path.join(root,'vendor/three.min.js'));
const source=fs.readFileSync(path.join(root,'viewer.js'),'utf8').replace('  function disposeGroup','  window.__auditLine=line;\n  function disposeGroup');
const context={window:{THREE:T},console};vm.createContext(context);vm.runInContext(source,context);
const make=context.window.__auditLine,front=make([[0,0,0],[0,.1,0]],0x7946b0),rear=make([[0,0,0],[0,.1,0]],0x7946b0,true,true);
if(front.material.depthWrite||rear.material.depthWrite)throw Error('An edge can overwrite face depth');
if(rear.material.depthFunc!==T.GreaterDepth||front.material.depthFunc!==T.LessEqualDepth)throw Error('Front/rear depth tests changed');
if(rear.material.opacity>=front.material.opacity)throw Error('Rear stroke is not lighter');
const baseRender=Object.getPrototypeOf(context.window.SectorView.prototype).render;
const scene=new T.Scene();scene.add(front,rear);
const camera=new T.OrthographicCamera(-2,2,2,-2,.01,100),rows=[];
const view={renderer:{render(){}},scene,camera,height:520,layoutLabels(){},renderMarkers(){}};
for(const zoom of [1,10,80]){
 camera.zoom=zoom;baseRender.call(view);
 const pxPerWorld=520*zoom/4,dash=rear.material.dashSize*pxPerWorld,gap=rear.material.gapSize*pxPerWorld;
 if(Math.abs(dash-6)>1e-12||Math.abs(gap-5)>1e-12)throw Error('Dashes expand with zoom');
 rows.push({zoom,nominalDashPixels:dash,nominalGapPixels:gap});
}
const result={status:'PASS',faceDepthPreserved:true,frontAndRearDepthTestsCorrect:true,rearOpacity:rear.material.opacity,frontOpacity:front.material.opacity,zoomChecks:rows};
fs.writeFileSync(path.join(__dirname,'line-depth-dash-regression.json'),JSON.stringify(result,null,2)+'\n');
process.stdout.write(JSON.stringify(result,null,2)+'\n');
