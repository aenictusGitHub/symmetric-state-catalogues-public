/* N=6 exact-graph projection, state inversion and nearest-face picking. */
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const root=path.resolve(__dirname,'../sector-length-sets');
const M=require(path.join(root,'majorana-math.js')),T=require(path.join(root,'vendor/three.min.js'));
const Occluder=require(path.join(root,'static-occlusion.js')),katex=require(path.join(root,'vendor/katex/katex.min.js'));
for(const file of ['n6-surface-bounds','n6-cyclic-constructors','n6-filtered-paired-constructors','n6-roof-constructors','n6-surface-state'])require(path.join(root,file+'.js'));
const G=JSON.parse(fs.readFileSync(path.join(root,'data/n6-smooth-geometry.json'))),full=process.argv.includes('--full'),randomMode=process.argv.includes('--random');
let seed=61207,count=0,worst=0,stars=0,adjusted=0,maxAdjustment=0;const random=()=>{seed=(1664525*seed+1013904223)>>>0;return seed/2**32;};
const close=(a,b,tol=1e-7)=>assert.ok(Math.abs(a-b)<=tol,`${a} != ${b}`),families=new Map(),fail=[];
for(const patch of G.patches){
  const ids=patch.attained.flatMap((yes,i)=>yes?[i]:[]),chosen=full?ids:ids.filter((_,i)=>i%Math.max(1,Math.floor(ids.length/100))===0);
  for(const i of chosen){
    const vertices=patch.triangles[i].map(j=>patch.positions[j]);
    const weights=full&&!randomMode?[1/3,1/3,1/3]:(()=>{const r=Math.sqrt(random()),s=random();return [1-r,r*(1-s),r*s];})();
    const raw=[0,1,2].map(k=>vertices.reduce((a,p,j)=>a+weights[j]*p[k],0));
    const state=N6SurfaceState.construct(raw,patch.surface,patch.kind,[0,1,2].map(k=>vertices.reduce((sum,p)=>sum+p[k],0)/3));
    if(!state){fail.push({kind:patch.kind,surface:patch.surface,point:N6SurfaceBounds.project(raw,patch.surface)});continue;}
    count++;worst=Math.max(worst,state.reconstructionResidual);if(state.surfaceAdjusted){adjusted++;maxAdjustment=Math.max(maxAdjustment,...state.surfaceAdjustment.map(Math.abs));const restored=N6SurfaceState.fromId(state.id);assert.ok(restored,'Adjusted point restores: '+state.id);restored.sectors.slice(0,3).forEach((v,j)=>close(v,state.sectors[j],1e-8));}
    close(state.coefficients.reduce((a,[r,i])=>a+r*r+i*i,0),1,2e-12);
    if(!families.has(state.family)){
      families.set(state.family,state);const e=state.analyticalExpression;
      for(const source of [state.mathLabel.slice(1,-1),e?.latex,e?.parameterLatex])if(source)katex.renderToString(source,{throwOnError:true});
      const result=M.majoranaStars(state.coefficients);assert.ok(result.reliable);assert.equal(result.stars.length,6);assert.ok(result.fidelity>1-1e-9);stars++;
      const restored=N6SurfaceState.fromId(state.id);assert.ok(restored);restored.sectors.slice(0,3).forEach((v,j)=>close(v,state.sectors[j],1e-8));
    }
  }
}
assert.equal(fail.length,0,JSON.stringify(fail.slice(0,20)));
// Audit graph formulas against the original double-precision export.
let graphError=0,graphs=0;
for(const patch of G.patches)if(['upper','lower'].includes(patch.surface))for(let i=0;i<patch.positions.length;i+=17){const[x,y,z]=patch.positions[i],got=N6SurfaceBounds[patch.surface](x,y);graphError=Math.max(graphError,Math.abs(got-z));close(got,z,1e-9);graphs++;}
const map=p=>[2*(p[0]+1)/9-1,2*(p[2]+1)/23-1,2*(p[1]+1.5)/17.5-1];
const surface=new T.Group();
for(const patch of G.patches){const values=[];for(const f of patch.triangles)for(const i of f)values.push(...map(patch.positions[i]));const geometry=new T.BufferGeometry();geometry.setAttribute('position',new T.Float32BufferAttribute(values,3));const mesh=new T.Mesh(geometry,new T.MeshBasicMaterial({side:T.DoubleSide}));mesh.userData.patch=patch;surface.add(mesh);}
surface.updateMatrixWorld(true);
const bvh=new Occluder(surface.children.map(m=>m.geometry.getAttribute('position').array));
const context={window:{THREE:T}};vm.runInNewContext(fs.readFileSync(path.join(root,'viewer.js'),'utf8'),context);
const View=context.window.SectorView,view=Object.create(View.prototype);
Object.assign(view,{n:6,data:{},width:600,height:600,surface,occlusionIndex:bvh,raycaster:new T.Raycaster(),camera:new T.OrthographicCamera(-1.8,1.8,1.8,-1.8,.01,100)});
let violet=0,gray=0,miss=0,rays=0;
for(const[yaw,pitch]of[[-1.04,.22],[-1.95,.50],[2.10,.377]]){
  view.camera.position.set(7*Math.cos(pitch)*Math.sin(yaw),7*Math.sin(pitch),7*Math.cos(pitch)*Math.cos(yaw));view.camera.lookAt(0,0,0);view.camera.updateMatrixWorld();
  for(let i=0;i<60;i++){
    const x=80+440*random(),y=80+440*random();view.raycaster.setFromCamera(new T.Vector2(2*x/600-1,1-2*y/600),view.camera);
    const ray=view.raycaster.ray,closest=bvh.closestHit(ray),brute=view.raycaster.intersectObjects(surface.children,false)[0];rays++;
    assert.equal(!!closest,!!brute);if(!closest){miss++;continue;}close(closest.distance,brute.distance,1e-7);
    assert.equal(bvh.occludes(ray,closest.distance+1e-5),true);assert.equal(bvh.occludes(ray,Math.max(closest.distance-1e-5,0)),false);
    const expected=G.patches[closest.source].attained[closest.faceIndex],hit=view.attainedSurfaceAt(x,y);
    assert.equal(!!hit,expected);if(expected){violet++;assert.equal(hit.surface,G.patches[closest.source].surface);const state=N6SurfaceState.construct(hit.point,hit.surface,hit.kind,hit.anchorPoint);assert.ok(state,'Visible violet hit reconstructs');}else gray++;
  }
}
assert.ok(violet>0&&gray>0);surface.visible=false;assert.equal(view.attainedSurfaceAt(300,300),null);
for(const id of ['n6-point:upper:NaN:6:8','n6-point:x0:0:5:100','n6-point:bad:0:5:0'])assert.equal(N6SurfaceState.fromId(id),null);
console.log(JSON.stringify({status:'PASS',full,randomMode,states:count,adjusted,maxAdjustment,families:[...families.keys()],worst,stars,graphs,graphError,rays,violet,gray,miss},null,2));
