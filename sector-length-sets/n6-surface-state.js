/* Verified inverses for the violet skin. Color alone never supplies a state. */
(function (g) {
  'use strict';
  function atPoint(point, surface, kind = '') {
    if (!g.N6SurfaceBounds || !g.N6SurfaceConstructors || !g.SectorMath) return null;
    const target = g.N6SurfaceBounds.project(point, surface);
    if (!target) return null;
    const [x,y,z] = target, tolerance = 1e-7;
    if (z < -tolerance || z > 20+tolerance || z < g.N6SurfaceBounds.lower(x,y)-tolerance || z > g.N6SurfaceBounds.upper(x,y)+tolerance) return null;
    const api=g.N6SurfaceConstructors;
    for (const name of ['filteredPlane','pairedFace','cyclic','roof']) {
      const state=api[name]?.(target,tolerance);
      if (!state) continue;
      const coefficients=g.SectorMath.canonicalState(state.coefficients);
      const actual=g.SectorMath.sectorLengths(coefficients);
      const residual=Math.max(...target.map((v,i)=>Math.abs(v-actual[i])));
      if (!Number.isFinite(residual) || residual>tolerance) continue;
      const stored=target;
      return {...state,coefficients,N:6,status:'attained',interactivePoint:true,
        label:'Six-qubit surface state',mathLabel:'$|\\psi(x,y,z)\\rangle$',
        id:`n6-point:${surface}:${stored.join(':')}`,
        sectors:actual,
        targetPoint:target,surface,surfaceKind:kind,reconstructionResidual:residual,
        proofOrigin:state.proofOrigin||state.source||state.family||'Explicit pure-state construction on an attained six-qubit surface',
        onOuterBoundary:true};
    }
    return null;
  }
  function construct(point, surface, kind = '', anchor = null) {
    const direct=atPoint(point,surface,kind);
    if(direct||!anchor)return direct;
    // Finite color triangles can cross the curved attainment seam. Retain
    // an actual reconstructed state throughout; never interpolate amplitudes
    // or interpret a failed inverse as an exclusion certificate.
    let result=atPoint(anchor,surface,kind);
    if(!result)return null;
    let lo=0,hi=1;
    for(let i=0;i<32;i++){
      const t=(lo+hi)/2,candidate=atPoint(point.map((v,k)=>anchor[k]+t*(v-anchor[k])),surface,kind);
      if(candidate){lo=t;result=candidate;}else hi=t;
    }
    // Stay just inside the last verified interval, avoiding tolerance-edge
    // sensitivity when the selected point is restored from a shared link.
    const interior=atPoint(point.map((v,k)=>anchor[k]+lo*(1-1e-5)*(v-anchor[k])),surface,kind);
    if(interior)result=interior;
    const clickedPoint=g.N6SurfaceBounds.project(point,surface);
    return {...result,surfaceAdjusted:true,clickedPoint,
      proofOrigin:'Nearby verified attained point at the surface color seam. '+result.proofOrigin,
      surfaceAdjustment:result.targetPoint.map((v,i)=>v-clickedPoint[i])};
  }
  function fromId(id){
    const match=/^n6-point:(upper|lower|x0|y5):([^:]+):([^:]+):([^:]+)$/.exec(id||'');
    return match?construct(match.slice(2).map(Number),match[1]):null;
  }
  g.N6SurfaceState=Object.freeze({construct,fromId});
})(globalThis);
