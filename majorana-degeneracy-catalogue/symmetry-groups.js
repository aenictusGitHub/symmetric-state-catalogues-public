/* Classify numerically detected proper SO(3) rotation groups.
 * Consumes MajoranaSymmetry.detect(...); never changes scientific state data. */
(function (root) {
  'use strict';
  let cache=new WeakMap();
  const disclaimer='Numerical classification of proper rotations only; reflections and inversion are excluded. This is not an exact symmetry proof.';

  function signature(axes) {
    const result={};
    for(const axis of axes) if(axis) result[axis.order]=(result[axis.order]||0)+1;
    return result;
  }

  function sameSignature(actual,expected) {
    const a=Object.keys(actual).sort(),b=Object.keys(expected).sort();
    return a.length===b.length && a.every((key,i)=>key===b[i] && actual[key]===expected[key]);
  }

  function gcd(a,b) {while(b) {[a,b]=[b,a%b];} return a;}

  function permutationOrder(p) {
    const seen=Array(p.length).fill(false);
    let order=1;
    for(let i=0;i<p.length;i++) {
      if(seen[i]) continue;
      let j=i,count=0;
      do {seen[j]=true;j=p[j];count++;} while(!seen[j]);
      order=order/gcd(order,count)*count;
    }
    return order;
  }

  function failure(kind,message,result) {
    return {status:kind==='invalid'?'invalid':'unclassified',kind,
      symbol:kind==='invalid'?'\\text{unavailable}':'\\text{unclassified}',
      label:kind==='invalid'?'Unavailable':'Unclassified',
      name:kind==='invalid'?'Symmetry unavailable':'Unclassified rotational symmetry',
      order:Number.isInteger(result?.rotationCount)?result.rotationCount:null,
      degree:null,continuous:result?.status==='continuous',numerical:true,
      properRotationsOnly:true,closureVerified:false,axisSignature:signature(Array.isArray(result?.axes)?result.axes:[]),
      tolerance:result?.tolerance??null,notes:[message,disclaimer]};
  }

  function finiteClosure(result) {
    const rotations=result.rotations,N=result.rotationCount,m=result.pointCount;
    if(!Array.isArray(rotations) || rotations.length!==N || !Number.isInteger(m) || m<2)
      return {passed:false,reason:'A complete finite rotation list and point count are required.'};
    const keys=new Set();
    for(const rotation of rotations) {
      const p=rotation.permutation;
      if(!Array.isArray(p) || p.length!==m || new Set(p).size!==m ||
         p.some(j=>!Number.isInteger(j)||j<0||j>=m))
        return {passed:false,reason:'The rotation list contains an invalid point permutation.'};
      const key=p.join(',');
      if(keys.has(key)) return {passed:false,reason:'The rotation list contains duplicate point permutations.'};
      if(rotation.order!==permutationOrder(p))
        return {passed:false,reason:'A rotation order disagrees with its point-permutation order.'};
      keys.add(key);
    }
    if(!keys.has(Array.from({length:m},(_,i)=>i).join(',')))
      return {passed:false,reason:'The identity rotation is missing.'};
    for(const a of rotations) for(const b of rotations) {
      const composition=b.permutation.map(j=>a.permutation[j]);
      if(!keys.has(composition.join(',')))
        return {passed:false,reason:'The detected permutations are not closed under composition at this tolerance.'};
    }
    // A finite set of permutations with identity and closure is a group, so
    // inverses follow without a separate numerical matrix-inverse tolerance.
    return {passed:true,reason:'The full finite permutation action contains identity and is closed under composition.'};
  }

  function classifyUncached(result) {
    if(!result || result.status==='invalid' || !Array.isArray(result.axes))
      return failure('invalid','No valid constellation symmetry result is available.',result);
    const axes=result.axes,axisSignature=signature(axes);
    let elementOrderSignature=null;
    const finish=(kind,symbol,label,name,degree,order,continuous,closureVerified,extraNotes=[])=>({
      status:'classified',kind,symbol,label,name,degree,order,continuous,numerical:true,
      properRotationsOnly:true,closureVerified,axisSignature,elementOrderSignature,tolerance:result.tolerance??null,
      notes:[...extraNotes,disclaimer]});

    if(result.status==='continuous') {
      if(axes.length!==1 || axes[0].order!==null || result.rotationCount!==null)
        return failure('unknown','The continuous-axis description is incomplete or inconsistent.',result);
      const families=result.continuousFamilies||[];
      if(families.length===0)
        return finish('axial-continuous','\\mathrm{SO}(2)','SO(2)',
          'Continuous axial rotation group',null,null,true,null,
          ['All rotation angles about one axis are allowed.']);
      if(families.length===1 && families[0].kind==='perpendicular-half-turn-axes' && families[0].order===2)
        return finish('dihedral-continuous','D_{\\infty}','D∞',
          'Continuous dihedral rotational group',null,null,true,null,
          ['The proper-rotation subgroup is SO(2) semidirect C2: arbitrary axial rotations and half-turns about every perpendicular axis.']);
      return failure('unknown','This continuous symmetry family is not one of the supported collinear cases.',result);
    }

    if(result.status!=='finite' || !Number.isInteger(result.rotationCount) || result.rotationCount<1)
      return failure('unknown','The finite rotation count is missing or invalid.',result);
    const N=result.rotationCount;
    if(axes.some(axis=>!axis||!Number.isInteger(axis.order)||axis.order<2||axis.continuous))
      return failure('unknown','Finite axes must have integral cyclic orders of at least2.',result);
    if(axes.some(axis=>axis.rotationCount!==axis.order-1) ||
       1+axes.reduce((sum,axis)=>sum+axis.order-1,0)!==N)
      return failure('unknown','The axis orders do not account for every nonidentity rotation exactly once.',result);
    const closure=finiteClosure(result);
    if(!closure.passed) return failure('unknown',closure.reason,result);
    elementOrderSignature={};
    for(const rotation of result.rotations)
      elementOrderSignature[rotation.order]=(elementOrderSignature[rotation.order]||0)+1;
    const impliedOrders={1:1};
    for(const axis of axes) for(let power=1;power<axis.order;power++) {
      const order=axis.order/gcd(axis.order,power);
      impliedOrders[order]=(impliedOrders[order]||0)+1;
    }
    if(!sameSignature(elementOrderSignature,impliedOrders))
      return {...failure('unknown','The element orders do not agree with the cyclic subgroups described by the axes.',result),closureVerified:true};
    const confirmed=[closure.reason];
    if(N===1 && axes.length===0)
      return finish('cyclic','C_{1}','C1','Trivial rotational group',1,1,false,true,confirmed);
    if(axes.length===1 && axes[0].order===N)
      return finish('cyclic',`C_{${N}}`,`C${N}`,'Cyclic rotational group',N,N,false,true,confirmed);
    if(N===12 && sameSignature(axisSignature,{2:3,3:4}))
      return finish('tetrahedral','T','T','Tetrahedral rotational group',null,12,false,true,confirmed);
    if(N===24 && sameSignature(axisSignature,{2:6,3:4,4:3}))
      return finish('octahedral','O','O','Octahedral rotational group',null,24,false,true,confirmed);
    if(N===60 && sameSignature(axisSignature,{2:15,3:10,5:6}))
      return finish('icosahedral','I','I','Icosahedral rotational group',null,60,false,true,confirmed);
    if(N%2===0) {
      const n=N/2,expected=n===2?{2:3}:{2:n,[n]:1};
      if(n>=2 && sameSignature(axisSignature,expected))
        return finish('dihedral',`D_{${n}}`,`D${n}`,'Dihedral rotational group',n,N,false,true,
          [...confirmed,`The notation D${n} denotes a proper rotation group of order ${2*n}.`]);
    }
    return {...failure('unknown','The closed rotation set does not match a supported SO(3) axis-order signature.',result),closureVerified:true};
  }

  function classify(result) {
    if(!result || typeof result!=='object') return classifyUncached(result);
    const key=JSON.stringify([result.status,result.rotationCount,result.pointCount,result.tolerance,
      result.axes?.map(a=>[a.order,a.rotationCount,a.continuous]),result.continuousFamilies,
      result.rotations?.map(r=>[r.order,r.permutation])]);
    const previous=cache.get(result);
    if(previous?.key===key) return previous.value;
    const value=classifyUncached(result);
    cache.set(result,{key,value});
    return value;
  }

  const api={classify,clearCache(){cache=new WeakMap();},version:'1.0.0'};
  root.MajoranaGroups=api;
  if(typeof module==='object' && module.exports) module.exports=api;
})(typeof window==='object'?window:globalThis);
