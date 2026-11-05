/* Numerical proper-rotation symmetries of a weighted Majorana constellation.
 * Offline, dependency-free. This is geometric verification at a tolerance,
 * not an exact symmetry or anticoherence proof. */
(function (root) {
  'use strict';

  const TAU = 2 * Math.PI;
  let cache = new WeakMap();
  const dot = (a, b) => a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  const cross = (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
  const length = a => Math.hypot(a[0], a[1], a[2]);
  const scale = (a, s) => a.map(x => x*s);
  const subtract = (a, b) => a.map((x, i) => x-b[i]);
  const unit = a => scale(a, 1/length(a));
  const distanceSquared = (a, b) => dot(subtract(a,b), subtract(a,b));
  const apply = (R, p) => [dot(R[0],p), dot(R[1],p), dot(R[2],p)];

  function canonicalAxis(a) {
    const v = unit(a);
    let pivot = 0;
    for (let k=1; k<3; k++) if (Math.abs(v[k]) > Math.abs(v[pivot])) pivot=k;
    return v[pivot] < 0 ? scale(v,-1) : v;
  }

  function frame(p, q) {
    const e1=p, e2=unit(subtract(q,scale(p,dot(p,q))));
    return [e1,e2,cross(e1,e2)]; // Orthonormal frame columns.
  }

  function frameRotation(source, target) {
    return [0,1,2].map(i => [0,1,2].map(j =>
      target[0][i]*source[0][j] + target[1][i]*source[1][j] + target[2][i]*source[2][j]));
  }

  function axisAngle(R) {
    // A diagonal quaternion branch avoids division by sin(theta) near pi.
    const tr=R[0][0]+R[1][1]+R[2][2];
    let w,x,y,z,s;
    if (tr > 0) {
      s=2*Math.sqrt(Math.max(0,1+tr));
      w=s/4; x=(R[2][1]-R[1][2])/s; y=(R[0][2]-R[2][0])/s; z=(R[1][0]-R[0][1])/s;
    } else if (R[0][0] >= R[1][1] && R[0][0] >= R[2][2]) {
      s=2*Math.sqrt(Math.max(0,1+R[0][0]-R[1][1]-R[2][2]));
      w=(R[2][1]-R[1][2])/s; x=s/4; y=(R[0][1]+R[1][0])/s; z=(R[0][2]+R[2][0])/s;
    } else if (R[1][1] >= R[2][2]) {
      s=2*Math.sqrt(Math.max(0,1+R[1][1]-R[0][0]-R[2][2]));
      w=(R[0][2]-R[2][0])/s; x=(R[0][1]+R[1][0])/s; y=s/4; z=(R[1][2]+R[2][1])/s;
    } else {
      s=2*Math.sqrt(Math.max(0,1+R[2][2]-R[0][0]-R[1][1]));
      w=(R[1][0]-R[0][1])/s; x=(R[0][2]+R[2][0])/s; y=(R[1][2]+R[2][1])/s; z=s/4;
    }
    const norm=Math.hypot(w,x,y,z);
    w/=norm; x/=norm; y/=norm; z/=norm;
    const v=Math.hypot(x,y,z), angle=2*Math.atan2(v,Math.abs(w));
    return {angle, direction:v>1e-13 ? canonicalAxis([x,y,z]) : null};
  }

  function bijection(R, points, tolerance) {
    const limit=tolerance*tolerance;
    const candidates=points.map(p => {
      const moved=apply(R,p.direction);
      return points.map((q,j) => ({j,d:distanceSquared(moved,q.direction),weight:q.weight}))
        .filter(q => q.weight===p.weight && q.d<=limit).sort((a,b)=>a.d-b.d);
    });
    if (candidates.some(row=>row.length===0)) return null;
    const targetOwner=Array(points.length).fill(-1), permutation=Array(points.length).fill(-1);
    function augment(i, seen) {
      for (const candidate of candidates[i]) {
        const j=candidate.j;
        if (seen[j]) continue;
        seen[j]=true;
        if (targetOwner[j]===-1 || augment(targetOwner[j],seen)) {
          targetOwner[j]=i; permutation[i]=j; return true;
        }
      }
      return false;
    }
    // Usually every row has one candidate. Augmentation also handles close stars.
    const order=points.map((_,i)=>i).sort((a,b)=>candidates[a].length-candidates[b].length);
    for (const i of order) if (!augment(i,Array(points.length).fill(false))) return null;
    const maxResidual=Math.sqrt(Math.max(...permutation.map((j,i)=>
      distanceSquared(apply(R,points[i].direction),points[j].direction))));
    return {permutation,maxResidual};
  }

  function permutationOrder(permutation) {
    const visited=Array(permutation.length).fill(false), cycles=[];
    for (let i=0; i<permutation.length; i++) {
      if (visited[i]) continue;
      let j=i, count=0;
      do {visited[j]=true; j=permutation[j]; count++;} while (!visited[j]);
      if (count>1) cycles.push(count);
    }
    // Every off-axis orbit of a nonidentity SO(3) rotation has the same length.
    if (new Set(cycles).size>1) return null;
    return cycles.length ? cycles[0] : 1;
  }

  function invalid(tolerance, message) {
    return {status:'invalid', axes:[], rotations:[], rotationCount:0,
      continuousFamilies:[], tolerance, maxResidual:null, notes:[message],
      scientificStatus:'Numerical geometric verification; not an exact symmetry proof.'};
  }

  function prepare(input,tolerance) {
    const raw=[];
    for (const item of input) {
      const a=Array.isArray(item) ? item.slice(0,3) : [item?.x,item?.y,item?.z];
      const weight=Array.isArray(item) ? (item[3] ?? 1) : (item?.multiplicity ?? item?.weight ?? 1);
      if (a.length!==3 || !a.every(Number.isFinite) || length(a)<1e-14 ||
          !Number.isInteger(weight) || weight<1) throw new Error('Points need finite nonzero coordinates and positive integer multiplicities.');
      raw.push({direction:unit(a),weight});
    }
    // Repeated rows and explicit weights are equivalent representations. Use a
    // much smaller merge tolerance than the symmetry matching tolerance.
    const mergeTolerance=Math.min(1e-12,tolerance/1000), points=[];
    for (const p of raw) {
      const q=points.find(q=>distanceSquared(p.direction,q.direction)<=mergeTolerance*mergeTolerance);
      if (q) {
        q.direction=unit(q.direction.map((x,k)=>x*q.weight+p.direction[k]*p.weight));
        q.weight+=p.weight;
      } else points.push({direction:p.direction.slice(),weight:p.weight});
    }
    let minimumSeparation=Infinity;
    for (let i=0;i<points.length;i++) for(let j=0;j<i;j++)
      minimumSeparation=Math.min(minimumSeparation,Math.sqrt(distanceSquared(points[i].direction,points[j].direction)));
    return {points,minimumSeparation,mergedInputPoints:input.length-points.length};
  }

  function compute(input,tolerance) {
    if (!Array.isArray(input) || input.length===0) return invalid(tolerance,'No constellation points are available.');
    if (!Number.isFinite(tolerance) || tolerance<=0 || tolerance>=0.1)
      return invalid(tolerance,'The matching tolerance must be positive and smaller than0.1.');
    let prepared;
    try {prepared=prepare(input,tolerance);} catch (error) {return invalid(tolerance,error.message);}
    const {points,minimumSeparation,mergedInputPoints}=prepared, count=points.length;
    const notes=['Proper SO(3) rotations only. Multiplicities are preserved. Results are numerical at the stated unit-sphere chord-distance tolerance.'];
    if (minimumSeparation<4*tolerance) notes.push('Some distinct points are close to the matching tolerance; the numerical symmetry assignment may be ambiguous.');
    if (mergedInputPoints) notes.push(`${mergedInputPoints} coincident input rows were combined into weighted points.`);
    const base={tolerance, pointCount:count, inputPointCount:input.length,
      totalMultiplicity:points.reduce((sum,p)=>sum+p.weight,0), minimumSeparation:Number.isFinite(minimumSeparation)?minimumSeparation:null,
      mergedInputPoints, notes, scientificStatus:'Numerical geometric verification; not an exact symmetry proof.'};

    let reference=null, bestCross=0;
    for(let i=0;i<count;i++) for(let j=i+1;j<count;j++) {
      const sine=length(cross(points[i].direction,points[j].direction));
      if(sine>bestCross) {bestCross=sine; reference=[i,j];}
    }
    if (count===1 || bestCross<=Math.max(1e-12,tolerance/4)) {
      const direction=canonicalAxis(points[0].direction), continuousFamilies=[];
      let positiveWeight=0,negativeWeight=0,maxResidual=0;
      for (const p of points) {
        if(dot(p.direction,direction)>=0) positiveWeight+=p.weight; else negativeWeight+=p.weight;
        // A full axial turn can move an almost-collinear point by twice its
        // perpendicular distance. This is an upper bound for every angle.
        maxResidual=Math.max(maxResidual,2*length(cross(p.direction,direction)));
      }
      if(positiveWeight===negativeWeight && positiveWeight>0) {
        continuousFamilies.push({kind:'perpendicular-half-turn-axes',normal:direction.slice(),
          order:2,rotationAngle:Math.PI,
          description:'Every axis perpendicular to the main axis gives a half-turn symmetry.'});
        notes.push('Equal antipodal weights have a continuous family of perpendicular half-turn axes, in addition to rotations about the main axis.');
      } else notes.push('Arbitrary rotation angles about the displayed axis preserve this collinear constellation.');
      return {...base,status:'continuous',rotationCount:null,rotations:[],continuousFamilies,maxResidual,
        axes:[{direction,order:null,continuous:true,angle:null,generatorAngle:null,rotationCount:null,maxResidual}]};
    }

    const [i,j]=reference, sourceFrame=frame(points[i].direction,points[j].direction);
    const gram=dot(points[i].direction,points[j].direction), rotations=[], seen=new Set();
    let testedCandidatePairs=0;
    for(let a=0;a<count;a++) {
      if(points[a].weight!==points[i].weight) continue;
      for(let b=0;b<count;b++) {
        if(a===b || points[b].weight!==points[j].weight ||
           Math.abs(dot(points[a].direction,points[b].direction)-gram)>4*tolerance) continue;
        if(length(cross(points[a].direction,points[b].direction))<bestCross/2) continue;
        testedCandidatePairs++;
        const R=frameRotation(sourceFrame,frame(points[a].direction,points[b].direction));
        const matched=bijection(R,points,tolerance);
        if(!matched) continue;
        const key=matched.permutation.join(',');
        if(seen.has(key)) continue;
        const order=permutationOrder(matched.permutation);
        if(order===null) continue;
        const extracted=axisAngle(R);
        if(order>1 && (!extracted.direction ||
          Math.abs(Math.sin(order*extracted.angle/2))>Math.max(1e-10,8*order*tolerance/bestCross))) continue;
        seen.add(key);
        rotations.push({matrix:R,angle:order===1?0:extracted.angle,order,
          axis:order===1?null:extracted.direction,permutation:matched.permutation,
          maxResidual:matched.maxResidual});
      }
    }
    rotations.sort((a,b)=>a.order-b.order || a.angle-b.angle);
    if(bestCross<Math.sqrt(tolerance)) notes.push('The constellation is nearly collinear; extracted finite axes are poorly conditioned.');
    const axes=[],axisTolerance=Math.min(1e-3,Math.max(1e-8,20*tolerance/bestCross));
    for(const rotation of rotations) {
      if(rotation.order===1) continue;
      let axis=axes.find(axis=>Math.min(distanceSquared(axis.direction,rotation.axis),
        distanceSquared(axis.direction,scale(rotation.axis,-1)))<=axisTolerance*axisTolerance);
      if(!axis) {
        axis={direction:rotation.axis.slice(),order:rotation.order,continuous:false,
          angle:TAU/rotation.order,generatorAngle:TAU/rotation.order,rotationCount:0,maxResidual:0,detectedAngles:[]};
        axes.push(axis);
      }
      axis.order=Math.max(axis.order,rotation.order); axis.angle=axis.generatorAngle=TAU/axis.order;
      axis.rotationCount++; axis.maxResidual=Math.max(axis.maxResidual,rotation.maxResidual);
      axis.detectedAngles.push(rotation.angle);
    }
    axes.sort((a,b)=>b.order-a.order || a.direction[0]-b.direction[0] || a.direction[1]-b.direction[1] || a.direction[2]-b.direction[2]);
    if(rotations.length===1) notes.push('No nonidentity proper rotation was detected at this tolerance.');
    return {...base,status:'finite',axes,rotations,rotationCount:rotations.length,
      continuousFamilies:[],maxResidual:Math.max(0,...rotations.map(r=>r.maxResidual)),
      referencePair:reference,referencePairSine:bestCross,testedCandidatePairs,
      method:'Oriented frame maps of a maximally nonparallel reference pair, followed by a complete multiplicity-preserving bijection.'};
  }

  function detect(points,options={}) {
    const tolerance=options.tolerance ?? 1e-6;
    if(!Array.isArray(points)) return compute(points,tolerance);
    const signature=JSON.stringify([tolerance,points]);
    const previous=cache.get(points);
    if(previous?.signature===signature) return previous.result;
    const result=compute(points,tolerance);
    cache.set(points,{signature,result});
    return result;
  }

  const api={detect,clearCache(){cache=new WeakMap();},version:'1.0.0'};
  root.MajoranaSymmetry=api;
  if(typeof module==='object' && module.exports) module.exports=api;
})(typeof window==='object'?window:globalThis);
