'use strict';
// Self-contained geometry/marker audit; executes production addPatch without a DOM or WebGL.
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert/strict');
const root=path.resolve(__dirname,'..');
const T=require(path.join(root,'vendor/three.min.js'));
const load=name=>JSON.parse(fs.readFileSync(path.join(root,'data',name),'utf8'));
const smooth=load('n6-smooth-geometry.json'),original=load('n6.json'),nodes=load('n6-manuscript-nodes.json');
const context={window:{THREE:T},AtlasMath:{state:(id,label)=>label,axis:text=>text},console};
context.window.AtlasMath=context.AtlasMath;
vm.runInNewContext(fs.readFileSync(path.join(root,'viewer.js'),'utf8'),context,{filename:'viewer.js'});
const view=Object.create(context.window.SectorView.prototype);
view.n=6;view.surface=new T.Group();
const mapped=p=>[2*(p[0]+1)/9-1,2*(p[2]+1)/23-1,2*(p[1]+1.5)/17.5-1];
const cross=(a,b,c)=>{const u=b.map((v,i)=>v-a[i]),v=c.map((v,i)=>v-a[i]);return [u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];};
const length=a=>Math.hypot(...a),dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const report={passed:true,patches:[],triangles:0,expandedVertexCount:0,sharedEdgeWindingErrors:0,nonmanifoldEdges:0,warnings:[]};
for(const [patchIndex,patch] of smooth.patches.entries()){
  assert.equal(patch.worldNormals.length,patch.positions.length,'one indexed world normal per vertex');
  assert.equal(patch.attained.length,patch.triangles.length,'one attainment classification per face');
  view.addPatch(patch);
  const geometry=view.surface.children.at(-1).geometry;
  const normal=geometry.getAttribute('normal').array,position=geometry.getAttribute('position').array;
  assert.equal(position.length,patch.triangles.length*9);
  assert.equal(normal.length,position.length);
  const positions=patch.positions.map(mapped),sums=positions.map(()=>[0,0,0]),faceNormals=[],faceAreas=[],edges=new Map();
  let outputOffset=0,maxLengthError=0;
  for(const [ti,tri] of patch.triangles.entries()){
    const n=cross(...tri.map(id=>positions[id])),area=length(n)/2;
    assert.ok(area>0&&Number.isFinite(area),'nondegenerate finite triangle');
    faceNormals.push(n.map(v=>v/(2*area)));faceAreas.push(area);
    for(const id of tri){
      assert.ok(id>=0&&id<positions.length,'valid triangle index');
      for(let axis=0;axis<3;axis++,outputOffset++){
        assert.equal(position[outputOffset],Math.fround(positions[id][axis]),'production expanded world position');
        assert.equal(normal[outputOffset],Math.fround(patch.worldNormals[id][axis]),'production expands indexed normal for the correct vertex');
        sums[id][axis]+=n[axis];
      }
    }
    for(let i=0;i<3;i++){
      const a=tri[i],b=tri[(i+1)%3],key=Math.min(a,b)+':'+Math.max(a,b);
      if(!edges.has(key))edges.set(key,[]);
      edges.get(key).push({triangle:ti,forward:a<b});
    }
  }
  let maxWeightedNormalError=0;
  for(let i=0;i<sums.length;i++){
    const len=length(sums[i]);
    const expected=len?sums[i].map(v=>v/len):[0,0,0];
    maxWeightedNormalError=Math.max(maxWeightedNormalError,...expected.map((v,j)=>Math.abs(v-patch.worldNormals[i][j])));
    maxLengthError=Math.max(maxLengthError,Math.abs(length(patch.worldNormals[i])-1));
  }
  assert.ok(maxLengthError<2e-8,'unit stored normals');
  assert.ok(maxWeightedNormalError<2e-8,'stored normals equal area-weighted indexed face normals');
  let windingErrors=0,nonmanifold=0,minAdjacentFaceDot=1;
  for(const fan of edges.values()){
    if(fan.length>2)nonmanifold++;
    if(fan.length===2){
      if(fan[0].forward===fan[1].forward)windingErrors++;
      minAdjacentFaceDot=Math.min(minAdjacentFaceDot,dot(faceNormals[fan[0].triangle],faceNormals[fan[1].triangle]));
    }
  }
  assert.equal(windingErrors,0,'consistent orientation along shared edges');
  assert.equal(nonmanifold,0,'no nonmanifold shared edges');
  const opposed=[];
  patch.triangles.forEach((tri,ti)=>{
    const minDot=Math.min(...tri.map(id=>dot(faceNormals[ti],patch.worldNormals[id])));
    if(minDot<-.001)opposed.push({triangle:ti,worldArea:faceAreas[ti],minimumNormalDot:minDot,sectorVertices:tri.map(id=>patch.positions[id])});
  });
  const worldArea=faceAreas.reduce((s,v)=>s+v,0),opposedArea=opposed.reduce((s,v)=>s+v.worldArea,0);
  const info={patchIndex,kind:patch.kind,surface:patch.surface,triangles:patch.triangles.length,positions:positions.length,
    maxStoredNormalLengthError:maxLengthError,maxAreaWeightedNormalError:maxWeightedNormalError,sharedEdgeWindingErrors:windingErrors,nonmanifoldEdges:nonmanifold,
    minimumAdjacentFaceDot:minAdjacentFaceDot,worldArea,opposedTriangleCount:opposed.length,opposedTriangleArea:opposedArea,opposedAreaFraction:opposedArea/worldArea};
  report.patches.push(info);report.triangles+=patch.triangles.length;report.expandedVertexCount+=patch.triangles.length*3;
  report.sharedEdgeWindingErrors+=windingErrors;report.nonmanifoldEdges+=nonmanifold;
  if(opposed.length)report.warnings.push({...info,opposedTriangles:opposed});
}
assert.equal(report.triangles,smooth.metadata.triangles);

// Run the actual marker grouping path with DOM/render-only operations stubbed.
const markerView=Object.create(context.window.SectorView.prototype);
markerView.target=new T.Vector3();markerView.content=new T.Group();
for(const name of ['surface','samples','edges','nodes'])markerView[name]=new T.Group();
markerView.clearLabels=function(){this.markers=[];this.labels=[];};
markerView.drawAxes=()=>{};markerView.resize=()=>{};markerView.applyOptions=()=>{};
markerView.addLabel=()=>({});
const named=original.states.filter(s=>nodes.some(n=>n.attained&&n.stateIds.includes(s.id)));
markerView.setData({...original,states:named,plotNodes:nodes,outer:{patches:[],edges:[],nodes}});
const violet=markerView.markers.filter(m=>m.kind==='sector-violet').length;
const gray=markerView.markers.filter(m=>m.kind==='sector-gray').length;
assert.equal(nodes.length,27);assert.equal(violet,16);assert.equal(gray,11);
assert.equal(markerView.markers.length,27,'production marker grouping preserves the original 27 manuscript nodes');
for(const n of nodes)assert.equal(markerView.markers.filter(m=>m.world.every((v,i)=>Math.abs(v-mapped(n.point)[i])<1e-8)).length,1,'exactly one marker at every manuscript tuple');
report.markers={originalManuscriptNodes:nodes.length,renderedMarkers:markerView.markers.length,violet,gray,allOriginalTuplesRepresentedExactlyOnce:true};
report.warningInterpretation='Opposed fan normals are flagged, not suppressed. Shared-edge index winding is consistent; some thin clipped geometric triangles bend sharply relative to their area-weighted vertex fan. Large patch creases are intentionally kept separate.';
fs.writeFileSync(path.join(__dirname,'smooth-normals-marker-regression.json'),JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify({passed:report.passed,triangles:report.triangles,expandedVertexCount:report.expandedVertexCount,markers:report.markers,sharedEdgeWindingErrors:report.sharedEdgeWindingErrors,nonmanifoldEdges:report.nonmanifoldEdges,warnings:report.warnings.map(({patchIndex,kind,surface,opposedTriangleCount,opposedTriangleArea,opposedAreaFraction})=>({patchIndex,kind,surface,opposedTriangleCount,opposedTriangleArea,opposedAreaFraction}))},null,2));
