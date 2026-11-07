/* Integration checks against the actual viewer renderer, without a browser.
 * Test-only VM instrumentation exposes its scene, never changes the public API.
 */
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const near=(a,b,eps=1e-10)=>assert.ok(Math.abs(a-b)<=eps*Math.max(1,Math.abs(a),Math.abs(b)),`${a} != ${b}`);
const escape=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function element(){
 const handlers=new Map(),classes=new Set();
 return {style:{},children:[],checked:false,value:'',textContent:'',innerHTML:'',
  handlers,classList:{add:c=>classes.add(c),remove:c=>classes.delete(c),contains:c=>classes.has(c),toggle:(c,on)=>on?classes.add(c):classes.delete(c)},
  addEventListener:(n,f)=>handlers.set(n,f),removeEventListener:(n,f)=>{if(handlers.get(n)===f)handlers.delete(n);},
  setAttribute(){},append(...items){this.children.push(...items);},replaceChildren(...items){this.children=[...items];},
  getBoundingClientRect(){return {width:486,height:parseFloat(this.style.height)||437};}
 };
}
const controls=new Map(),get=id=>{if(!controls.has(id))controls.set(id,element());return controls.get(id);};
const help=element(),frame=element();frame.querySelector=()=>help;
const gradient={addColorStop(){}};
const context=new Proxy({measureText:s=>({width:s.length*7}),createRadialGradient:()=>gradient},{get:(o,k)=>k in o?o[k]:()=>{}});
const canvas=get('canvas');canvas.parentElement=frame;canvas.getContext=()=>context;
get('fit').checked=true;get('shape').value='sphere';get('order').value='highest';
const window={
 SubspaceSphereLayout:require('../subspace-sphere-layout.js'),
 HusimiRendering:require('../husimi-rendering.js'),
 SubspaceRenderGeometry:require('../subspace-render-geometry.js'),
 SubspaceSvgExport:{mathLabel:(tex,x,y)=>`<text x="${x}" y="${y}" data-tex="${escape(tex)}">${escape(tex)}</text>`},
 katex:{renderToString:tex=>`<span>${escape(tex)}</span>`},devicePixelRatio:1
};
const document={createElement:tag=>{const e=element();if(tag==='canvas')e.getContext=()=>context;return e;}};
const sourcePath=path.join(__dirname,'../subspace-viewer.js');let source=fs.readFileSync(sourcePath,'utf8');
const hook='return {align,set(next,r,c,nextMesh)';assert.equal(source.split(hook).length,2,'Renderer API changed; review the test hook.');
source=source.replace(hook,'return {sceneForTest:scene,renderForTest:render,align,set(next,r,c,nextMesh)');
const publicHook='window.SubspaceViewer={open,close};';assert.equal(source.split(publicHook).length,2);
source=source.replace(publicHook,publicHook+'window.RendererForTest=makeRenderer;window.DialogForTest=d=>dialog=d;');
let frameId=0;
vm.runInNewContext(source,{window,document,requestAnimationFrame:()=>++frameId,cancelAnimationFrame(){},ResizeObserver:class{observe(){}disconnect(){}},setTimeout,Map,console},{filename:sourcePath});
window.DialogForTest({querySelector:selector=>get(selector.replace('#sv-',''))});
const records=[{rank:2,weight:1},{rank:4,weight:.5},{rank:6,weight:.1}].map(r=>({...r,available:true,selectedPoints:[[1,0,0]],oppositePoints:[[-1,0,0]]}));
const model={n:6,multipoles:records.map(r=>({rank:r.rank,weight:r.weight,zero:false}))};
const mesh={qMin:0,qMax:1,vertices:[{n:[1,0,0],q:0},{n:[0,1,0],q:.2},{n:[0,0,1],q:.5},{n:[-1,0,0],q:1}],faces:[[0,1,2,3]],contours:{segments:[]}};
const renderer=window.RendererForTest(model,mesh,{axes:[]});
const scene=()=>renderer.sceneForTest(486,437,true);
const shells=objects=>objects.filter(o=>o.type==='shell');
const radii=objects=>shells(objects).map(o=>o.r);
renderer.set('operator','all',{all:true,records});
const nested=scene();assert.equal(nested.filter(o=>o.type==='clip-start').length,0);assert.equal(shells(nested).length,3);
assert.equal(new Set(shells(nested).map(o=>`${o.cx},${o.cy}`)).size,1,'Nested shells must share a center');
near(radii(nested)[1]/radii(nested)[0],.5);near(radii(nested)[2]/radii(nested)[0],.1);
get('separate').checked=true;
const separated=scene(),r=radii(separated);assert.equal(separated.filter(o=>o.type==='clip-start').length,3);assert.equal(separated.filter(o=>o.type==='clip-end').length,3);
assert.equal(shells(separated).length,3);assert.equal(new Set(shells(separated).map(o=>`${o.cx},${o.cy}`)).size,3,'Separate shells need separate centers');
near(r[1]/r[0],.5);near(r[2]/r[0],.1);
// Fit-selected remains checked: separation must still preserve the norm radii.
assert.ok(r[0]>r[1]&&r[1]>r[2],'Separate spheres must ignore per-rank fit enlargement');
for(const rank of [2,4,6]){
 assert.ok(separated.some(o=>o.tex===`\\ell=${rank}`),'Per-sphere rank label absent');
 assert.ok(separated.some(o=>o.tex?.startsWith(`w_{${rank}}`)),'Per-sphere norm label absent');
}
let activeClip=null;
for(const o of separated){
 if(o.type==='clip-start')activeClip=o;
 else if(o.type==='clip-end')activeClip=null;
 else if(o.type==='shell'){
  assert.ok(activeClip,'Shell outside card clipping group');
  assert.ok(o.cx-o.r>=activeClip.x-1e-9&&o.cx+o.r<=activeClip.x+activeClip.width+1e-9);
  assert.ok(o.cy-o.r>=activeClip.y-1e-9&&o.cy+o.r<=activeClip.y+activeClip.height+1e-9);
 }
}
const svg=renderer.exportSvg();
const clipIds=[...svg.matchAll(/<clipPath id="([^"]+)"/g)].map(m=>m[1]);
assert.equal(clipIds.length,3);assert.equal(new Set(clipIds).size,3);
assert.equal((svg.match(/<g clip-path=/g)||[]).length,3);
assert.equal((svg.match(/<radialGradient /g)||[]).length,3);
assert.ok(!/NaN|Infinity/.test(svg));
for(const id of clipIds)assert.ok(svg.includes(`clip-path="url(#${id})"`));
// Scroll layout is intrinsic canvas height; leaving it restores the normal view.
renderer.renderForTest();assert.equal(frame.classList.contains('sv-separated'),true);assert.ok(parseFloat(canvas.style.height)>=310);
assert.equal(get('math-labels').style.height,canvas.style.height);
assert.equal(get('math-labels').children.filter(c=>c.className==='sv-math-cell').length,3,'Mathematical overlays need the same cell clipping groups');
for(const c of get('math-labels').children)assert.equal(c.children.length,2,'Each card needs its own rank and norm labels');
const wheel=canvas.handlers.get('wheel');let prevented=false;
wheel({deltaY:40,preventDefault(){prevented=true;}});assert.equal(prevented,false,'Plain scrolling must pass to the gallery');
wheel({deltaY:0,ctrlKey:true,preventDefault(){prevented=true;}});assert.equal(prevented,true,'Ctrl+scroll should retain zoom');
get('separate').checked=false;renderer.renderForTest();assert.equal(canvas.style.height,'');assert.equal(frame.classList.contains('sv-separated'),false);
const restored=scene();assert.deepEqual(radii(restored),radii(nested));
renderer.set('operator',4,records[1]);near(radii(scene())[0],radii(nested)[0]);
get('fit').checked=false;near(radii(scene())[0],radii(nested)[0]*.5);
get('separate').checked=true;renderer.set('husimi',4,null);
const husimi=scene();assert.equal(husimi.filter(o=>o.type==='clip-start').length,0);assert.equal(husimi.filter(o=>o.type==='poly').length,1);
renderer.dispose();
const manyRecords=Array.from({length:13},(_,i)=>({...records[0],rank:i+1,weight:(i+1)/13}));
const many=window.RendererForTest({n:13,multipoles:manyRecords.map(r=>({rank:r.rank,weight:r.weight,zero:false}))},mesh,{axes:[]});
get('separate').checked=true;get('coordinate').checked=true;
many.set('operator','all',{all:true,records:manyRecords});many.renderForTest();
assert.equal(canvas.style.height,'950px','Thirteen ranks should create a genuinely scrollable gallery');
assert.equal(get('math-labels').style.height,'950px');
const overlays=get('math-labels').children;assert.equal(overlays.length,13);
for(const c of overlays){
 assert.equal(c.className,'sv-math-cell');assert.equal(c.children.length,5,'Coordinates, rank and weight remain inside the card overlay');
 const width=Number(c.style.cssText.match(/width:([\d.]+)px/)[1]);
 for(const label of c.children){const left=Number(label.style.cssText.match(/left:([\d.-]+)px/)[1]);assert.ok(left>=0&&left<width,'Math label origin escaped its card');}
}
many.set('operator',7,manyRecords[6]);many.renderForTest();assert.equal(canvas.style.height,'');assert.equal(frame.classList.contains('sv-separated'),false);
many.dispose();

// A signed field must use a neutral zero, and negative radial lobes must
// remain in their original directions when magnitude display is selected.
const dirs=[[0,0,-1],[.1,0,-Math.sqrt(.99)],[0,.1,-Math.sqrt(.99)],[-.1,0,-Math.sqrt(.99)]];
const signedMesh={...mesh,qMin:-.2,qMax:1,colorMax:1,vertices:dirs.map(n=>({n,q:-.2})),faces:[[0,1,2,3]],contours:{segments:[]}};
const positiveMesh={...signedMesh,qMin:.2,qMax:1,vertices:dirs.map(n=>({n,q:.2}))};
const signed=window.RendererForTest(model,signedMesh,{axes:[{axis:[0,0,1],order:4}]});
get('separate').checked=false;get('coordinate').checked=false;get('shape').value='radial';
signed.set('wigner',2,null,signedMesh);
const negativeGeometry=signed.sceneForTest(486,437,true).filter(o=>o.type==='poly');
signed.set('husimi',2,null,positiveMesh);
const positiveGeometry=signed.sceneForTest(486,437,true).filter(o=>o.type==='poly');
assert.deepEqual(negativeGeometry.map(o=>o.points),positiveGeometry.map(o=>o.points),'Negative Wigner magnitudes must not flip to the antipode');
get('shape').value='sphere';get('axes').checked=true;
const patches={...signedMesh,qMin:-1,vertices:[...dirs.map(n=>({n,q:-1})),...dirs.map(n=>({n,q:0})),...dirs.map(n=>({n,q:1}))],faces:[[0,1,2,3],[4,5,6,7],[8,9,10,11]]};
signed.set('wigner',2,null,patches);
const patchColors=signed.sceneForTest(486,437,true).filter(o=>o.type==='poly').map(o=>o.fill);
assert.ok(patchColors.includes('#fbfaf4'),'Zero Wigner value must be neutral white');
assert.ok(patchColors.includes('#2b5895'),'Negative Wigner value must be blue');
assert.ok(patchColors.includes('#b1411c'),'Positive Wigner value must be orange');
const signedSvg=signed.exportSvg();assert.ok(signedSvg.includes('data-tex="C_{4}"'),'Cn guide must be LaTeX in exports');
assert.ok(!signedSvg.includes('<image'),'Wigner export must remain vector geometry');
signed.dispose();

// Refined sphere contours form one continuous visible path. The rear half of
// a closed level loop must be removed rather than connected across the sphere.
const level=.5,loop=[[-.8,.6,0],[0,1,0],[.8,.6,0],[.8,-.6,0],[0,-1,0],[-.8,-.6,0]];
const loopSegments=loop.map((a,i)=>({a,b:loop[(i+1)%loop.length],level}));
const radialA=[-.15,Math.sqrt(1-.15**2-.1**2),.1],radialB=[.15,Math.sqrt(1-.15**2-.1**2),.1];
const cap=[[-.5,Math.SQRT1_2,-.5],[.5,Math.SQRT1_2,-.5],[.5,Math.SQRT1_2,.5],[-.5,Math.SQRT1_2,.5]];
const refinedMesh={qMin:-1,qMax:1,colorMax:1,vertices:cap.map(n=>({n,q:level})),faces:[[0,1,2,3]],
 contours:{levels:[level],segments:loopSegments},radialContours:{levels:[level],segments:[{a:radialA,b:radialB,level}]}};
get('contours').checked=true;get('axes').checked=false;get('shape').value='sphere';
const continuous=window.RendererForTest(model,refinedMesh,{axes:[]});continuous.set('wigner',2,null,refinedMesh);
const contourPaths=continuous.sceneForTest(486,437,true).filter(o=>o.type==='path');
assert.equal(contourPaths.length,1,'The visible half-loop should be one joined path');
assert.ok(contourPaths[0].points.length>=5,'Joining must preserve the visible contour vertices');
assert.ok(contourPaths[0].points.every(p=>p.z>=-1e-12),'Sphere contours must never reconnect through the hidden hemisphere');
near(contourPaths[0].points[0].z,0);near(contourPaths[0].points.at(-1).z,0);
assert.equal(contourPaths[0].width,1);assert.equal(contourPaths[0].opacity,1);
assert.ok(!['#fff','#ffffff','#faf8f1'].includes(contourPaths[0].color),'Contours must use a single stroke without a white halo');
const continuousSvg=continuous.exportSvg();
assert.equal((continuousSvg.match(/stroke-linejoin="round"/g)||[]).length,1,'SVG and canvas must contain the same joined contour path');
assert.ok(!continuousSvg.includes('<image'),'Joined contours must remain vector paths in SVG');
// Radial clipping uses the contour geometry matched to the rendered coarse
// surface, not the separately refined colored-sphere geometry.
get('shape').value='radial';
const radialPaths=continuous.sceneForTest(486,437,true).filter(o=>o.type==='path');
assert.equal(radialPaths.length,1);assert.equal(radialPaths[0].points.length,2);
function projectDefault(v){
 const [x,y,z]=v,cx=Math.cos(-.45),sx=Math.sin(-.45),cy=Math.cos(-.25),sy=Math.sin(-.25),rx=cx*x+sx*z,rz=-sx*x+cx*z,scale=Math.min(486*.39,437*.40);
 return {x:486/2+scale*rx,y:437/2-scale*(sy*y+cy*rz),z:cy*y-sy*rz};
}
const expectedRadial=[radialA,radialB].map(v=>projectDefault(v.map(x=>x*level)));
for(const expected of expectedRadial)assert.ok(radialPaths[0].points.some(actual=>Math.hypot(actual.x-expected.x,actual.y-expected.y,actual.z-expected.z)<1e-9),'Radial mode must select the original surface-matched contour');
const radialSvg=continuous.exportSvg();assert.equal((radialSvg.match(/stroke-linejoin="round"/g)||[]).length,radialPaths.length);
assert.ok(!/NaN|Infinity/.test(radialSvg));continuous.dispose();

// Scalar axes retain their front/back style, but portions hidden by the
// opaque body are faint. Check the actual rendered sections and SVG opacities.
get('contours').checked=false;get('axes').checked=true;get('coordinate').checked=false;get('shape').value='sphere';get('order').value='highest';
const sphereAxes=window.RendererForTest(model,{...mesh,colorMax:1},{axes:[{axis:[0,0,1],order:4}]});
const axisLines=objects=>objects.filter(o=>o.type==='line'&&typeof o.axisOccluded==='boolean');
const cssScale=Math.min(486*.39,437*.40);
const physicalRadius=p=>Math.hypot((p.x-486/2)/cssScale,(p.y-437/2)/cssScale,p.z);
function checkOpacitySections(lines){
 assert.ok(lines.some(o=>o.axisOccluded),'An axis must visibly pass inside/behind the body');
 assert.ok(lines.some(o=>!o.axisOccluded),'The exposed axis ends must retain their original opacity');
 for(const line of lines){
  assert.ok(line.dash.length===0||line.dash.join(',')==='5,5','Rear dashes must be preserved');
  const base=line.dash.length?.45:.95;
  near(line.opacity,base*(line.axisOccluded?.22:1));
 }
}
sphereAxes.set('husimi',2,null);
const sphereSections=axisLines(sphereAxes.sceneForTest(486,437,true));checkOpacitySections(sphereSections);
const frontExposed=sphereSections.find(o=>!o.axisOccluded&&!o.dash.length),rearExposed=sphereSections.find(o=>!o.axisOccluded&&o.dash.length);
assert.ok(frontExposed&&rearExposed,'Both tilted axis ends extend beyond the sphere silhouette');
near(Math.min(physicalRadius(frontExposed.a),physicalRadius(frontExposed.b)),1,1e-9,'The front axis must emerge at the unit-sphere surface');
const defaultAxisDepth=projectDefault([0,0,1]).z,rearEmergence=1/Math.sqrt(1-defaultAxisDepth**2);
near(Math.min(physicalRadius(rearExposed.a),physicalRadius(rearExposed.b)),rearEmergence,1e-9);
assert.ok(rearEmergence>1,'A rear segment outside the physical sphere can still lie behind its silhouette');
assert.ok(rearExposed.dashOffset<0,'Splitting a rear guide must continue its original dash phase');
const axisSvg=sphereAxes.exportSvg(),svgAxisPaths=[...axisSvg.matchAll(/<path\b[^>]*stroke="#a45635"[^>]*>/g)].map(m=>m[0]);
assert.equal(svgAxisPaths.length,sphereSections.length,'Canvas and SVG must draw the same exposed/hidden sections');
for(const line of sphereSections){
 const matching=svgAxisPaths.find(path=>Math.abs(Number(path.match(/opacity="([^"]+)"/)[1])-line.opacity)<1e-12&&path.includes(`stroke-dasharray="${line.dash.join(' ')}"`));
 assert.ok(matching,'SVG must preserve each axis section opacity and rear dash style');
}
sphereAxes.set('wigner',2,null,{...mesh,colorMax:1});checkOpacitySections(axisLines(sphereAxes.sceneForTest(486,437,true)));
// Translucent operator shells keep their previous guides, not scalar-body
// occlusion. This also protects the separate clickable sphere gallery.
sphereAxes.set('operator',4,records[1]);
const operatorGuide=sphereAxes.sceneForTest(486,437,true).filter(o=>o.type==='line'&&o.color==='#a45635');
assert.ok(operatorGuide.length>0);assert.ok(operatorGuide.every(o=>o.axisOccluded!==true));for(const line of operatorGuide)near(line.opacity,line.dash.length?.45:.95);sphereAxes.dispose();

// A prolate radial surface has a much narrower waist than a unit sphere.
// Exposed portions must follow that shape even when contour lines are off.
const bodyVertices=[],bodyFaces=[],bodyRows=24,bodyCols=48;
for(let i=0;i<=bodyRows;i++)for(let j=0;j<=bodyCols;j++){
 const theta=Math.PI*i/bodyRows,phi=2*Math.PI*j/bodyCols,n=[Math.sin(theta)*Math.cos(phi),Math.sin(theta)*Math.sin(phi),Math.cos(theta)];
 bodyVertices.push({n,q:.4+.6*n[2]**2});
}
for(let i=0;i<bodyRows;i++)for(let j=0;j<bodyCols;j++){const a=i*(bodyCols+1)+j,b=a+bodyCols+1;bodyFaces.push([a,a+1,b+1,b]);}
const bodyMesh={qMin:.4,qMax:1,colorMax:1,vertices:bodyVertices,faces:bodyFaces,contours:{segments:[]}};
const bodyAxes=window.RendererForTest(model,bodyMesh,{axes:[{axis:[1,0,0],order:2}]});get('shape').value='radial';
bodyAxes.set('husimi',2,null,bodyMesh);
const radialSections=axisLines(bodyAxes.sceneForTest(486,437,true));checkOpacitySections(radialSections);
assert.ok(radialSections.filter(o=>!o.axisOccluded).some(o=>Math.min(physicalRadius(o.a),physicalRadius(o.b))<.8),'Radial guides must emerge at the actual waist rather than a unit-sphere cutoff');
assert.equal(get('contours').checked,false,'Axis occlusion must not depend on the contour checkbox');
const beforeRotation=radialSections.map(o=>[o.a.x,o.a.y,o.b.x,o.b.y]);
get('right').handlers.get('click')();
const afterRotation=axisLines(bodyAxes.sceneForTest(486,437,true));checkOpacitySections(afterRotation);
assert.notDeepEqual(afterRotation.map(o=>[o.a.x,o.a.y,o.b.x,o.b.y]),beforeRotation,'Axis visibility must be recomputed during camera rotation');
const radialAxisSvg=bodyAxes.exportSvg();assert.ok(!/NaN|Infinity/.test(radialAxisSvg));
assert.ok(radialAxisSvg.includes('opacity="0.209"')||radialAxisSvg.includes('opacity="0.099"'),'SVG must carry the faded radial sections');
get('axes').checked=false;assert.equal(axisLines(bodyAxes.sceneForTest(486,437,true)).length,0);bodyAxes.dispose();

console.log('PASS: actual nested/separate scenes, norm-radius ratios, per-card labels/SVG/DOM clipping, tall gallery scrolling, normal-view restoration and fit behavior, joined contour paths, horizon clipping, surface-matched radial contours and depth-aware scalar axis transparency.');
