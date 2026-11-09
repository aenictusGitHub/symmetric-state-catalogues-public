/*
 * Proper rotational symmetries of an unordered, multiplicity-weighted
 * Majorana constellation. Plain JavaScript; no graphics/network dependency.
 * Directions use SectorMath's Bloch coordinates, not the viewer's camera frame.
 *
 * For a non-axial constellation, a rotation is uniquely fixed by the images
 * of one ordered pair of non-collinear stars. Enumerating all admissible
 * image pairs and checking the ENTIRE weighted constellation therefore finds
 * every proper rotation, within the explicitly returned numerical tolerance.
 * Reflections and antipodal inversion are not counted as rotations.
 */
(function(root){
  'use strict';
  const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
  const scale=(a,s)=>a.map(v=>v*s);
  const sub=(a,b)=>a.map((v,i)=>v-b[i]);
  const distance=(a,b)=>Math.hypot(...sub(a,b));
  const unit=a=>scale(a,1/Math.hypot(...a));
  const clamp=v=>Math.max(-1,Math.min(1,v));
  const apply=(R,a)=>R.map(row=>dot(row,a));
  function canonicalAxis(a){
    a=unit(a);let k=0;for(let j=1;j<3;j++)if(Math.abs(a[j])>Math.abs(a[k]))k=j;
    return scale(a,a[k]<0?-1:1).map(v=>Math.abs(v)<1e-13?0:v);
  }
  function frame(a,b){
    const u=unit(a),v=unit(sub(b,scale(u,dot(u,b)))),w=cross(u,v);
    return [u,v,w]; // Orthonormal columns, with positive orientation.
  }
  function rotation(F,G){
    return [0,1,2].map(i=>[0,1,2].map(j=>G[0][i]*F[0][j]+G[1][i]*F[1][j]+G[2][i]*F[2][j]));
  }
  function sameMatrix(A,B,tolerance){
    return Math.hypot(...A.flatMap((r,i)=>r.map((v,j)=>v-B[i][j])))<tolerance;
  }
  function permutationError(R,points,multiplicities,tolerance){
    const used=new Set();let worst=0;
    for(let i=0;i<points.length;i++){
      const p=apply(R,points[i]);let best=-1,error=Infinity;
      for(let j=0;j<points.length;j++)if(!used.has(j)&&multiplicities[i]===multiplicities[j]){
        const e=distance(p,points[j]);if(e<error){best=j;error=e;}
      }
      if(best<0||error>tolerance)return null;
      used.add(best);worst=Math.max(worst,error);
    }
    return worst;
  }
  function axisAndAngle(R){
    const angle=Math.acos(clamp((R[0][0]+R[1][1]+R[2][2]-1)/2));
    if(angle<1e-7)return null;
    let axis;
    if(Math.abs(Math.sin(angle))>1e-5){
      axis=[R[2][1]-R[1][2],R[0][2]-R[2][0],R[1][0]-R[0][1]];
    }else{
      // At pi the skew part is zero. Recover aa^T=(R+I)/2 instead.
      let k=0;for(let j=1;j<3;j++)if(R[j][j]>R[k][k])k=j;
      axis=[0,0,0];axis[k]=Math.sqrt(Math.max(0,(R[k][k]+1)/2));
      for(let j=0;j<3;j++)if(j!==k)axis[j]=(R[j][k]+R[k][j])/(4*axis[k]);
    }
    return {direction:canonicalAxis(axis),angle};
  }
  function analyze(groups,options={}){
    const tolerance=options.tolerance??2e-6;
    if(!Array.isArray(groups)||!groups.length)throw new Error('A nonempty Majorana constellation is required.');
    const points=groups.map(g=>{
      const p=Array.isArray(g)?g:[g.x,g.y,g.z];
      if(p.length!==3||p.some(v=>!Number.isFinite(v))||Math.hypot(...p)<1e-12)throw new Error('Invalid Majorana direction.');
      return unit(p);
    });
    const multiplicities=groups.map(g=>g.multiplicity??1);
    if(multiplicities.some(m=>!Number.isInteger(m)||m<1))throw new Error('Star multiplicities must be positive integers.');
    const N=multiplicities.reduce((a,b)=>a+b,0);
    const common={tolerance,qubits:N,numerical:true};
    const antipodal=points.length===2&&Math.hypot(...points[0].map((v,j)=>v+points[1][j]))<tolerance;
    if(points.length===1||antipodal){
      const transverseC2=antipodal&&multiplicities[0]===multiplicities[1];
      return {...common,axes:[{direction:canonicalAxis(points[0]),order:'infinity'}],rotationCount:null,
        continuous:{axial:true,transverseC2},reliable:true,
        summary:transverseC2?'Continuous axial symmetry; also a continuous family of transverse C2 axes.':'Continuous axial symmetry.'};
    }
    // Choose the best-conditioned ordered source pair.
    let pair=null,strength=0;
    for(let i=0;i<points.length;i++)for(let j=i+1;j<points.length;j++){
      const s=Math.hypot(...cross(points[i],points[j]));if(s>strength){strength=s;pair=[i,j];}
    }
    if(!pair||strength<20*tolerance)return {...common,axes:[],rotationCount:1,continuous:null,reliable:false,
      summary:'The constellation is too close to a collinear limit to resolve its axes reliably.'};
    const [a,b]=pair,F=frame(points[a],points[b]),referenceDot=dot(points[a],points[b]);
    const rotations=[],axes=[];let maximumPermutationError=0;
    for(let i=0;i<points.length;i++)for(let j=0;j<points.length;j++){
      if(i===j||multiplicities[i]!==multiplicities[a]||multiplicities[j]!==multiplicities[b])continue;
      if(Math.abs(dot(points[i],points[j])-referenceDot)>2*tolerance)continue;
      const R=rotation(F,frame(points[i],points[j]));
      const error=permutationError(R,points,multiplicities,tolerance);
      if(error===null||rotations.some(S=>sameMatrix(R,S,10*tolerance)))continue;
      maximumPermutationError=Math.max(maximumPermutationError,error);rotations.push(R);
      const result=axisAndAngle(R);if(!result)continue;
      let order=0;
      // A nonidentity rotation of a non-axial N-star multiset has a cycle of
      // length equal to its order, hence order <= N. This also handles angles
      // 4pi/5, etc.; rounding 2pi/angle would give the wrong answer there.
      for(let k=2;k<=N;k++)if(Math.abs(Math.sin(k*result.angle/2))<20*k*tolerance){order=k;break;}
      if(!order)continue;
      const previous=axes.find(v=>Math.abs(dot(v.direction,result.direction))>1-1e-8);
      if(previous)previous.order=Math.max(previous.order,order);
      else axes.push({direction:result.direction,order});
    }
    axes.sort((a,b)=>b.order-a.order||b.direction[0]-a.direction[0]||b.direction[1]-a.direction[1]||b.direction[2]-a.direction[2]);
    return {...common,axes,rotationCount:rotations.length,continuous:null,reliable:true,maximumPermutationError,
      summary:axes.length?'Each displayed axis is verified to preserve every star and its multiplicity.':'No nontrivial rotational symmetry axis detected.'};
  }
  const api={analyze};root.MajoranaSymmetry=api;
  if(typeof module!=='undefined'&&module.exports)module.exports=api;
})(typeof globalThis!=='undefined'?globalThis:window);
