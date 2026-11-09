/*
 * Sector Atlas — dependency-free symmetric-state mathematics.
 *
 * Coefficients: c[k] = [real, imaginary] in the NORMALIZED Dicke basis
 * |psi> = sum(k=0..N) c[k] |D_N^(k)>, with k excitations |1>.
 *
 * Majorana convention: constituent |v> = a|0>+b|1>, z=b/a.
 * Roots of sum_k (-1)^k sqrt(binomial(N,k)) c[k] z^(N-k)
 * are the constituent spinors. Infinity is |1>; zero is |0>.
 * Bloch vectors are (2 Re z, 2 Im z, 1-|z|^2)/(1+|z|^2).
 * This is a constituent-state convention, not its antipodal dual.
 *
 * Polynomial roots are a visualization, not an algebraic certificate.
 * Every returned constellation is independently reconstructed and its
 * state fidelity is reported. Exact poles and coherent multiplicities
 * are handled without a numerical root solve.
 */
(function (root) {
  'use strict';
  const C = (r=0,i=0) => [r,i];
  const add=(a,b)=>[a[0]+b[0],a[1]+b[1]];
  const sub=(a,b)=>[a[0]-b[0],a[1]-b[1]];
  const mul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]];
  const conj=a=>[a[0],-a[1]];
  const scale=(a,s)=>[a[0]*s,a[1]*s];
  const abs2=a=>a[0]*a[0]+a[1]*a[1];
  const abs=a=>Math.hypot(a[0],a[1]);
  const div=(a,b)=>{const d=abs2(b);return [(a[0]*b[0]+a[1]*b[1])/d,(a[1]*b[0]-a[0]*b[1])/d];};
  function binomial(n,k) {if(k<0||k>n)return 0;let b=1;for(let i=1;i<=Math.min(k,n-k);i++)b=b*(n-i+1)/i;return b;}
  function asComplex(v) {return typeof v==='number'?[v,0]:[Number(v[0]),Number(v[1]||0)];}
  function normalizeState(coefficients) {
    if(!Array.isArray(coefficients)||coefficients.length<2) throw new Error('A state needs N+1 Dicke coefficients.');
    const c=coefficients.map(asComplex);
    if(c.some(z=>!Number.isFinite(z[0])||!Number.isFinite(z[1])))throw new Error('State coefficients must be finite.');
    const norm=Math.sqrt(c.reduce((s,z)=>s+abs2(z),0));
    if(!(norm>0))throw new Error('The zero vector is not a state.');
    return c.map(z=>scale(z,1/norm));
  }
  function canonicalState(coefficients, tolerance=1e-12) {
    const c=normalizeState(coefficients), first=c.find(z=>abs(z)>tolerance);
    const phase=scale(conj(first),1/abs(first));
    return c.map(z=>{const q=mul(z,phase);return [Math.abs(q[0])<tolerance?0:q[0],Math.abs(q[1])<tolerance?0:q[1]];});
  }
  function fidelity(a,b) {
    a=normalizeState(a); b=normalizeState(b);
    if(a.length!==b.length)return 0;
    let overlap=C();for(let k=0;k<a.length;k++)overlap=add(overlap,mul(conj(a[k]),b[k]));
    return Math.min(1,Math.max(0,abs2(overlap)));
  }
  function reducedDensity(coefficients,s) {
    const c=normalizeState(coefficients), n=c.length-1;
    if(!Number.isInteger(s)||s<0||s>n)throw new Error('Invalid reduction size.');
    const B=Array.from({length:s+1},(_,a)=>Array.from({length:n-s+1},(_,b)=>scale(c[a+b],Math.sqrt(binomial(s,a)*binomial(n-s,b)/binomial(n,a+b)))));
    return B.map((row,a)=>B.map((_,ap)=>{
      let v=C();for(let b=0;b<=n-s;b++)v=add(v,mul(row[b],conj(B[ap][b])));return v;
    }));
  }
  function reducedPurity(coefficients,s) {
    return reducedDensity(coefficients,s).reduce((sum,row)=>sum+row.reduce((v,z)=>v+abs2(z),0),0);
  }
  function sectorLengths(coefficients) {
    const c=normalizeState(coefficients), n=c.length-1, sectors=[1];
    for(let s=1;s<=n;s++){
      let v=2**s*reducedPurity(c,s);
      for(let k=0;k<s;k++)v-=binomial(s,k)*sectors[k]/binomial(n,k);
      v*=binomial(n,s);
      sectors.push(Math.abs(v)<2e-11?0:v);
    }
    return sectors.slice(1);
  }
  // Descending complex polynomial coefficients.
  function monic(p,tolerance=2e-13) {
    let a=p.map(asComplex), max=Math.max(...a.map(abs),0);
    if(!max)return [C()];
    while(a.length>1&&abs(a[0])<=tolerance*max)a.shift();
    if(a.length===1)return abs(a[0])<=tolerance*max?[C()]:[C(1)];
    const lead=a[0];return a.map(z=>div(z,lead));
  }
  function derivative(p){const n=p.length-1;return n?p.slice(0,-1).map((z,k)=>scale(z,n-k)):[C()];}
  function evaluate(p,z){let v=C();for(const a of p)v=add(mul(v,z),a);return v;}
  function polynomialDivision(a,b) {
    const r=a.map(z=>z.slice()), n=a.length-b.length;
    if(n<0)return {quotient:[C()],remainder:r};
    const q=Array.from({length:n+1},()=>C());
    for(let k=0;k<=n;k++){
      q[k]=div(r[k],b[0]);
      for(let j=0;j<b.length;j++)r[k+j]=sub(r[k+j],mul(q[k],b[j]));
    }
    return {quotient:q,remainder:r.slice(n+1)};
  }
  function gcd(a,b) {
    a=monic(a);b=monic(b);
    for(let step=0;step<12&&b.length>1;step++){
      const r=polynomialDivision(a,b).remainder;
      if(!r.length||Math.max(...r.map(abs))<2e-11)return b;
      a=b;b=monic(r,2e-11);
    }
    return [C(1)];
  }
  function squareFree(p) {
    let g=gcd(p,derivative(p)),w=monic(polynomialDivision(p,g).quotient),i=1;
    const factors=[];
    while(w.length>1&&i<=p.length){
      const y=gcd(w,g), f=monic(polynomialDivision(w,y).quotient);
      if(f.length>1)factors.push({polynomial:f,multiplicity:i});
      w=y;g=monic(polynomialDivision(g,y).quotient);i++;
    }
    if(factors.reduce((s,f)=>s+(f.polynomial.length-1)*f.multiplicity,0)!==p.length-1)return [{polynomial:p,multiplicity:1}];
    return factors;
  }
  function shifted(p,center){
    let out=[p[0]];
    for(let k=1;k<p.length;k++){
      const next=Array.from({length:out.length+1},()=>C());
      for(let j=0;j<out.length;j++){next[j]=add(next[j],out[j]);next[j+1]=add(next[j+1],mul(out[j],center));}
      next[next.length-1]=add(next[next.length-1],p[k]);out=next;
    }
    return out;
  }
  function simpleRoots(input) {
    const p=monic(input),n=p.length-1;
    if(n===0)return [];
    if(n===1)return [scale(p[1],-1)];
    const center=scale(p[1],-1/n), q=shifted(p,center);
    let radius=0;for(let k=1;k<=n;k++)radius=Math.max(radius,abs(q[k])**(1/k));
    if(radius<1e-14)return Array.from({length:n},()=>center.slice());
    const a=q.map((v,k)=>scale(v,radius**(-k))),da=derivative(a);
    let best=[],bestError=Infinity;
    for(let attempt=0;attempt<5;attempt++){
      let roots=Array.from({length:n},(_,i)=>{const angle=2*Math.PI*(i+.37+attempt*.137)/n,rr=1.15+attempt*.11;return C(rr*Math.cos(angle),rr*Math.sin(angle));});
      for(let iteration=0;iteration<250;iteration++){
        let maxStep=0;const next=[];
        for(let i=0;i<n;i++){
          let z=roots[i],v=evaluate(a,z),d=evaluate(da,z);
          if(abs(d)<1e-22){z=add(z,C(1e-9*(i+1),1e-9*(iteration+1)));v=evaluate(a,z);d=evaluate(da,z);}
          const newton=div(v,d);let sum=C();
          for(let j=0;j<n;j++)if(i!==j){let sep=sub(z,roots[j]);if(abs(sep)<1e-22)sep=C(1e-12*(i+1),1e-12*(j+1));sum=add(sum,div(C(1),sep));}
          let step=div(newton,sub(C(1),mul(newton,sum)));
          if(!Number.isFinite(step[0])||!Number.isFinite(step[1]))step=scale(newton,.25);
          if(abs(step)>3)step=scale(step,3/abs(step));
          next.push(sub(z,step));maxStep=Math.max(maxStep,abs(step));
        }
        roots=next;if(maxStep<3e-14)break;
      }
      roots=roots.map(z=>{for(let i=0;i<8;i++){const d=evaluate(da,z);if(abs(d)<1e-20)break;const step=div(evaluate(a,z),d);z=sub(z,step);if(abs(step)<1e-15)break;}return z;});
      const error=Math.max(...roots.map(z=>abs(evaluate(a,z))));
      if(error<bestError){bestError=error;best=roots;}
      if(error<1e-12)break;
    }
    return best.map(z=>add(center,scale(z,radius)));
  }
  function rootStar(z) {
    const r=Math.hypot(z[0],z[1]);
    if(r>1e8){const inv=1/r,den=1+inv*inv;return {x:2*(z[0]/r)*inv/den,y:2*(z[1]/r)*inv/den,z:(inv*inv-1)/den};}
    const den=1+r*r;return {x:2*z[0]/den,y:2*z[1]/den,z:(1-r*r)/den};
  }
  function starSpinor(star) {
    const len=Math.hypot(star.x,star.y,star.z);if(!(len>0))throw new Error('Invalid Majorana star.');
    const x=star.x/len,y=star.y/len,z=star.z/len;
    // Use the hemisphere with a well-conditioned nonzero spinor component.
    // sqrt(1+z) alone loses relative precision close to the south pole.
    let a,b;
    if(z>=0){a=C(Math.sqrt((1+z)/2));b=C(x/(2*a[0]),y/(2*a[0]));}
    else {b=C(Math.sqrt((1-z)/2));a=C(x/(2*b[0]),-y/(2*b[0]));}
    const norm=Math.sqrt(abs2(a)+abs2(b));return [scale(a,1/norm),scale(b,1/norm)];
  }
  function convolve(a,b){const c=Array.from({length:a.length+b.length-1},()=>C());for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)c[i+j]=add(c[i+j],mul(a[i],b[j]));return c;}
  function linearPower(a,b,n){let p=[C(1)];for(let k=0;k<n;k++)p=convolve(p,[a,b]);return p;}
  function rotatedState(c,a,b){
    const n=c.length-1,p=Array.from({length:n+1},()=>C());
    // U=[[a,-conj(b)],[b,conj(a)]]. The transformed Majorana form is
    // sum q_k (a*z-b)^(N-k) (conj(b)*z+conj(a))^k.
    for(let k=0;k<=n;k++){
      const term=convolve(linearPower(scale(b,-1),a,n-k),linearPower(conj(a),conj(b),k));
      const q=scale(c[k],(k%2?-1:1)*Math.sqrt(binomial(n,k)));
      for(let j=0;j<=n;j++)p[j]=add(p[j],mul(q,term[j]));
    }
    return normalizeState(Array.from({length:n+1},(_,k)=>scale(p[n-k],(k%2?-1:1)/Math.sqrt(binomial(n,k)))));
  }
  function undoStarRotation(star,a,b){
    const [ap,bp]=starSpinor(star),aa=add(mul(conj(a),ap),mul(conj(b),bp)),bb=sub(mul(a,bp),mul(b,ap));
    const overlap=mul(conj(aa),bb),norm=abs2(aa)+abs2(bb);
    return {x:2*overlap[0]/norm,y:2*overlap[1]/norm,z:(abs2(aa)-abs2(bb))/norm};
  }
  function stateFromStars(stars) {
    let p=[C(1)];
    for(const star of stars){
      const [a,b]=starSpinor(star);
      const next=Array.from({length:p.length+1},()=>C());
      for(let k=0;k<p.length;k++){next[k]=add(next[k],mul(p[k],a));next[k+1]=add(next[k+1],mul(p[k],b));}p=next;
    }
    const n=stars.length;return normalizeState(p.map((v,k)=>scale(v,1/Math.sqrt(binomial(n,k)))));
  }
  function groupedStars(stars,tolerance=2e-5) {
    const groups=[];
    for(const s of stars){const g=groups.find(g=>Math.hypot(g.x-s.x,g.y-s.y,g.z-s.z)<tolerance);
      if(g)g.multiplicity++;else groups.push({...s,multiplicity:1});}
    return groups;
  }
  function constellationResidual(c,stars){
    const n=c.length-1,polynomial=c.map((z,k)=>scale(z,(k%2?-1:1)*Math.sqrt(binomial(n,k))));
    let residual=0;
    for(const s of stars){
      const [a,b]=starSpinor(s);let v=C();
      for(let k=0;k<=n;k++){let term=polynomial[k];for(let j=0;j<n-k;j++)term=mul(term,b);for(let j=0;j<k;j++)term=mul(term,a);v=add(v,term);}
      residual=Math.max(residual,abs(v));
    }
    return residual;
  }
  function majoranaStars(coefficients,internalChart=false) {
    const c=normalizeState(coefficients),n=c.length-1,rho=reducedDensity(c,1);
    const bloch={x:2*rho[0][1][0],y:-2*rho[0][1][1],z:rho[0][0][0]-rho[1][1][0]};
    const dipole=Math.hypot(bloch.x,bloch.y,bloch.z);
    let stars=[],method='square-free Aberth roots';
    if(dipole>1-5e-13){
      const s={x:bloch.x/dipole,y:bloch.y/dipole,z:bloch.z/dipole};stars=Array.from({length:n},()=>({...s}));method='coherent multiplicity';
    }else{
      let p=c.map((z,k)=>scale(z,(k%2?-1:1)*Math.sqrt(binomial(n,k))));
      const tol=2e-14;let south=0,north=0;
      while(p.length>1&&abs(p[0])<tol){p.shift();south++;}
      while(p.length>1&&abs(p[p.length-1])<tol){p.pop();north++;}
      const polys=squareFree(monic(p)), roots=[];
      for(const f of polys)for(const r of simpleRoots(f.polynomial))for(let m=0;m<f.multiplicity;m++)roots.push(r.slice());
      stars=[...Array.from({length:north},()=>({x:0,y:0,z:1})),...roots.map(rootStar),...Array.from({length:south},()=>({x:0,y:0,z:-1}))];
      // Approximate square-free decomposition is accepted only after a separate
      // physical-state reconstruction. Otherwise retain all simple roots.
      if(stars.length!==n||fidelity(c,stateFromStars(stars))<1-2e-10){
        stars=[...Array.from({length:north},()=>({x:0,y:0,z:1})),...simpleRoots(monic(p)).map(rootStar),...Array.from({length:south},()=>({x:0,y:0,z:-1}))];method='Aberth roots';
      }
    }
    let reconstructed=stateFromStars(stars), f=fidelity(c,reconstructed),residual=constellationResidual(c,stars);
    // A tiny leading coefficient can place one root extremely far from all
    // others. Change the collective SU(2) chart, solve there, and rotate back.
    // This does not rotate the displayed physical constellation.
    if(!internalChart&&(f<1-2e-14||residual>1e-9)){
      const charts=[];
      for(let k=0;k<18;k++){
        const z=1-2*(k+.5)/18,theta=Math.acos(z),phi=k*Math.PI*(3-Math.sqrt(5));
        const a=C(Math.cos(theta/2)),b=C(Math.sin(theta/2)*Math.cos(phi),Math.sin(theta/2)*Math.sin(phi));
        const cc=rotatedState(c,a,b),score=Math.min(abs(cc[0]),abs(cc[n]));charts.push({a,b,cc,score});
      }
      charts.sort((a,b)=>b.score-a.score);
      for(const chart of charts.slice(0,4)){
        const trial=majoranaStars(chart.cc,true),candidate=trial.stars.map(s=>undoStarRotation(s,chart.a,chart.b));
        const candidateFidelity=fidelity(c,stateFromStars(candidate)),candidateResidual=constellationResidual(c,candidate);
        if(candidateFidelity>f+1e-15||(candidateFidelity>=f-2e-15&&candidateResidual<residual)){
          stars=candidate;f=candidateFidelity;residual=candidateResidual;method='adaptive SU(2) chart; '+trial.method;
        }
        if(f>1-2e-14&&residual<1e-12)break;
      }
    }
    stars.sort((a,b)=>b.z-a.z||Math.atan2(a.y,a.x)-Math.atan2(b.y,b.x));
    return {stars,groups:groupedStars(stars),fidelity:f,residual,method,reliable:f>1-1e-8};
  }
  const api={binomial,normalizeState,canonicalState,fidelity,reducedDensity,reducedPurity,sectorLengths,majoranaStars,stateFromStars,groupedStars,
    convention:'Normalized Dicke order k=0..N; constituent-state stars; |0> north, |1> south; z=b/a.'};
  root.SectorMath=api;
  if(typeof module!=='undefined'&&module.exports)module.exports=api;
})(typeof globalThis!=='undefined'?globalThis:window);
