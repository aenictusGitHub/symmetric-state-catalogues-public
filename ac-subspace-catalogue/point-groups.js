/* Numerical full O(3) point groups of weighted Majorana constellations.
 * Consumes MajoranaSymmetry and MajoranaGroups results. Reflections and
 * inversion are geometric symmetries here, not additional spin rotations. */
(function(root) {
  'use strict';
  let cache=new WeakMap();
  const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
  const norm=a=>Math.hypot(...a),scale=(a,s)=>a.map(x=>x*s);
  const unit=a=>scale(a,1/norm(a));
  const subtract=(a,b)=>a.map((x,i)=>x-b[i]);
  const distance=(a,b)=>norm(subtract(a,b));
  const apply=(M,p)=>M.map(row=>dot(row,p));
  const multiply=(A,B)=>A.map(row=>[0,1,2].map(j=>row.reduce((s,x,k)=>s+x*B[k][j],0)));
  const transpose=M=>[0,1,2].map(i=>M.map(row=>row[i]));
  const identity=[[1,0,0],[0,1,0],[0,0,1]],inversion=identity.map(row=>row.map(x=>-x));
  const matrixDistance=(A,B)=>Math.max(...A.flatMap((row,i)=>row.map((x,j)=>Math.abs(x-B[i][j]))));
  const disclaimer='Numerical classification of the full weighted geometric point group in O(3). Improper operations are not spin rotations. This is not an exact symmetry or anticoherence proof.';

  function prepare(input,tolerance) {
    const points=[],mergeTolerance=Math.min(1e-9,tolerance/10);
    for(const p of input) {
      const vector=Array.isArray(p)?p.slice(0,3):[p?.x,p?.y,p?.z];
      const weight=Array.isArray(p)?(p[3]??1):(p?.multiplicity??p?.weight??1);
      if(vector.length!==3 || !vector.every(Number.isFinite) || norm(vector)<1e-14 || !Number.isInteger(weight) || weight<1)
        throw new Error('Finite nonzero coordinates and positive integer multiplicities are required.');
      const direction=unit(vector),same=points.find(q=>distance(q.direction,direction)<=mergeTolerance);
      if(same){same.direction=unit(same.direction.map((x,i)=>x*same.weight+direction[i]*weight));same.weight+=weight;}
      else points.push({direction,weight});
    }
    return points;
  }

  function match(M,points,tolerance) {
    const moved=points.map(p=>apply(M,p.direction));
    const candidates=points.map((p,i)=>points.map((q,j)=>({j,d:distance(moved[i],q.direction),weight:q.weight}))
      .filter(q=>q.weight===p.weight&&q.d<=tolerance).sort((a,b)=>a.d-b.d));
    if(candidates.some(row=>!row.length))return null;
    const owners=Array(points.length).fill(-1),permutation=Array(points.length).fill(-1);
    function augment(i,visited) {
      for(const {j} of candidates[i]) {
        if(visited[j])continue;visited[j]=true;
        if(owners[j]===-1||augment(owners[j],visited)){owners[j]=i;permutation[i]=j;return true;}
      }
      return false;
    }
    for(const i of points.map((_,i)=>i).sort((a,b)=>candidates[a].length-candidates[b].length))
      if(!augment(i,Array(points.length).fill(false)))return null;
    return {permutation,maxResidual:Math.max(...permutation.map((j,i)=>distance(moved[i],points[j].direction)))};
  }

  function frame(p,q) {
    const e1=p,e2=unit(subtract(q,scale(p,dot(p,q))));
    return [e1,e2,cross(e1,e2)];
  }

  function improperMap(source,target) {
    // Map the selected pair and reverse its oriented normal: det(M)=-1.
    return [0,1,2].map(i=>[0,1,2].map(j=>
      target[0][i]*source[0][j]+target[1][i]*source[1][j]-target[2][i]*source[2][j]));
  }

  function findImproper(points,tolerance) {
    let reference=null,bestCross=0;
    for(let i=0;i<points.length;i++)for(let j=i+1;j<points.length;j++){
      const sine=norm(cross(points[i].direction,points[j].direction));
      if(sine>bestCross){bestCross=sine;reference=[i,j];}
    }
    if(!reference||bestCross<1e-12)return {matrix:null,bestCross,reference};
    const [i,j]=reference,source=frame(points[i].direction,points[j].direction),gram=dot(points[i].direction,points[j].direction);
    for(let a=0;a<points.length;a++){
      if(points[a].weight!==points[i].weight)continue;
      for(let b=0;b<points.length;b++){
        if(a===b||points[b].weight!==points[j].weight||
          Math.abs(dot(points[a].direction,points[b].direction)-gram)>4*tolerance||
          norm(cross(points[a].direction,points[b].direction))<bestCross/2)continue;
        const M=improperMap(source,frame(points[a].direction,points[b].direction)),verified=match(M,points,tolerance);
        if(verified)return {matrix:M,...verified,bestCross,reference};
      }
    }
    return {matrix:null,bestCross,reference};
  }

  function classifyUncached(input,rotations,rotationalGroup) {
    const group=rotationalGroup||root.MajoranaGroups?.classify(rotations);
    const tolerance=rotations?.tolerance??1e-6;
    const base={numerical:true,properRotationsOnly:false,rotationalSymbol:group?.symbol??null,
      rotationalOrder:group?.order??null,tolerance,notes:[disclaimer]};
    const failure=message=>({...base,status:'unclassified',kind:'unknown',symbol:'\\text{unclassified}',label:'Unclassified',
      name:'Full point group not determined',degree:null,order:null,continuous:false,includesImproper:null,
      improperOperationCount:null,maxResidual:null,closureVerified:false,notes:[message,disclaimer]});
    if(!Array.isArray(input)||!input.length||!rotations||group?.status!=='classified')
      return failure('A valid constellation and classified proper-rotation subgroup are required.');
    if(!Number.isFinite(tolerance)||tolerance<=0||tolerance>=.1)return failure('The matching tolerance is invalid.');
    let points;try{points=prepare(input,tolerance);}catch(error){return failure(error.message);}
    if(points.length!==rotations.pointCount)return failure('Constellation points do not match the proper-rotation result.');
    const finish=(kind,symbol,label,name,degree,order,continuous,includesImproper,extra={})=>({
      ...base,status:'classified',kind,symbol,label,name,degree,order,continuous,includesImproper,
      improperOperationCount:continuous?null:includesImproper?group.order:0,
      maxResidual:rotations.maxResidual??null,closureVerified:continuous?null:group.closureVerified,...extra});
    if(group.kind==='axial-continuous')
      return finish('axial-continuous','C_{\\infty v}','C∞v','Continuous axial point group',null,null,true,true,
        {notes:['All axial rotations and all reflection planes through the axis preserve this collinear constellation.',disclaimer]});
    if(group.kind==='dihedral-continuous')
      return finish('dihedral-continuous','D_{\\infty h}','D∞h','Continuous dihedral point group',null,null,true,true,
        {notes:['Equal antipodal weights also allow inversion and reflection in the perpendicular plane.',disclaimer]});
    if(rotations.status!=='finite'||!Array.isArray(rotations.rotations)||!rotations.rotations.every(r=>Array.isArray(r.matrix)))
      return failure('Complete proper-rotation matrices are required for the finite full point group.');
    const found=findImproper(points,tolerance);
    if(!found.matrix)return finish(group.kind,group.symbol,group.label,group.name.replace('rotational group','point group'),
      group.degree,group.order,false,false,{notes:['No orientation-reversing weighted symmetry was found; the full group equals its proper-rotation subgroup.',disclaimer]});

    // An improper symmetry generates the entire determinant-minus-one coset.
    // Verify every member against all points, not just the reference pair.
    const improper=rotations.rotations.map(r=>multiply(r.matrix,found.matrix)),matches=improper.map(M=>match(M,points,tolerance));
    if(matches.some(x=>!x))return failure('The improper coset does not verify consistently at this numerical tolerance.');
    const matrixTolerance=Math.max(1e-9,20*tolerance/found.bestCross);
    const isProper=M=>rotations.rotations.some(r=>matrixDistance(M,r.matrix)<=matrixTolerance);
    if(!isProper(multiply(found.matrix,found.matrix))||
      rotations.rotations.some(r=>!isProper(multiply(multiply(found.matrix,r.matrix),transpose(found.matrix)))))
      return failure('The orientation-reversing candidate does not generate a consistent group extension.');
    const extra={maxResidual:Math.max(rotations.maxResidual||0,...matches.map(x=>x.maxResidual)),closureVerified:true,
      improperWitness:found.matrix,referencePair:found.reference,
      method:'Orientation-reversing reference-pair frames; weighted bijections for the complete improper coset; matrix group-extension checks.'};
    const done=(kind,symbol,label,name,degree=null)=>finish(kind,symbol,label,name,degree,2*group.order,false,true,extra);
    const hasInversion=Boolean(match(inversion,points,tolerance));
    if(group.kind==='icosahedral')return done('icosahedral','I_h','Ih','Full icosahedral point group');
    if(group.kind==='octahedral')return done('octahedral','O_h','Oh','Full octahedral point group');
    if(group.kind==='tetrahedral')return hasInversion?
      done('tetrahedral','T_h','Th','Tetrahedral point group with inversion'):
      done('tetrahedral','T_d','Td','Full tetrahedral point group');
    const n=group.degree;
    if(group.kind==='cyclic'&&n===1)return hasInversion?
      done('inversion','C_i','Ci','Inversion point group'):
      done('reflection','C_s','Cs','Single-mirror point group');
    const principal=rotations.axes.filter(a=>a.order===n).map(a=>unit(a.direction));
    if(!principal.length)return failure('A principal rotation axis is missing.');
    const horizontal=principal.some(axis=>{
      const M=identity.map((row,i)=>row.map((x,j)=>x-2*axis[i]*axis[j]));
      return Boolean(match(M,points,tolerance));
    });
    if(group.kind==='dihedral')return horizontal?
      done('dihedral',`D_{${n}h}`,`D${n}h`,'Dihedral point group with a horizontal mirror',n):
      done('dihedral',`D_{${n}d}`,`D${n}d`,'Dihedral point group with diagonal mirrors',n);
    if(group.kind==='cyclic'){
      const vertical=improper.some(M=>distance(apply(M,principal[0]),principal[0])<=matrixTolerance);
      if(vertical)return done('cyclic',`C_{${n}v}`,`C${n}v`,'Cyclic point group with vertical mirrors',n);
      if(horizontal)return done('cyclic',`C_{${n}h}`,`C${n}h`,'Cyclic point group with a horizontal mirror',n);
      return done('improper-cyclic',`S_{${2*n}}`,`S${2*n}`,'Rotoreflection point group',2*n);
    }
    return failure('This proper-rotation subgroup is not supported for full-group classification.');
  }

  function classify(points,rotationsResult,rotationalGroup) {
    if(!Array.isArray(points))return classifyUncached(points,rotationsResult,rotationalGroup);
    const signature=JSON.stringify([points,rotationsResult,rotationalGroup]);
    const previous=cache.get(points);if(previous?.signature===signature)return previous.result;
    const result=classifyUncached(points,rotationsResult,rotationalGroup);cache.set(points,{signature,result});return result;
  }
  const api={classify,clearCache(){cache=new WeakMap();},version:'1.0.0'};
  root.MajoranaPointGroups=api;
  if(typeof module==='object'&&module.exports)module.exports=api;
})(typeof window==='object'?window:globalThis);
