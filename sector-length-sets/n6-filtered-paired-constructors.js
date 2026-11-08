/* Constructive N=6 face inverses, evaluated in floating point.
 * Sources: Julia/N6_filtered_octa_attainment.jl and
 * output/pdf/assets/N6_complete_paired_face_20261021.tex.
 * A successful inverse is independently checked with reduced-state purities.
 * Failure supplies no exclusion claim. No coefficient interpolation is used.
 */
(function(global){
  'use strict';
  const clamp=(x,a,b)=>Math.max(a,Math.min(b,x)), SQ=Math.sqrt;
  const zero=()=>Array.from({length:7},()=>[0,0]);
  const near=(a,b,tol)=>Math.max(...a.map((v,i)=>Math.abs(v-b[i])))<=tol;
  const number=x=>Number(x.toPrecision(10)).toString();
  function target(point,tolerance){return Array.isArray(point)&&point.length===3&&point.every(Number.isFinite)&&Number.isFinite(tolerance)&&tolerance>0;}
  function accepted(point,coefficients,family,parameters,analyticalExpression,tolerance){
    if(!coefficients||coefficients.length!==7||!coefficients.every(z=>z.length===2&&z.every(Number.isFinite)))return null;
    const norm=SQ(coefficients.reduce((s,z)=>s+z[0]*z[0]+z[1]*z[1],0));
    if(!(norm>0))return null;
    const c=coefficients.map(z=>[z[0]/norm,z[1]/norm]);
    const actual=global.SectorMath.sectorLengths(c);
    if(!near(actual.slice(0,3),point,tolerance))return null;
    return {coefficients:c,family,parameters,analyticalExpression};
  }
  const basisNote=String.raw`D_k=|D_6^{(k)}\rangle`;

  // Symmetric sixth power of the real single-qubit y-rotation.
  const choose=(n,k)=>k<0||k>n?0:[1,1,2,6,24,120,720][n]/([1,1,2,6,24,120,720][k]*[1,1,2,6,24,120,720][n-k]);
  function rotateY(c,beta){
    const a=Math.cos(beta/2),b=Math.sin(beta/2),out=Array(7).fill(0);
    for(let l=0;l<=6;l++)if(c[l]!==0)for(let k=0;k<=6;k++){
      let sum=0;
      for(let j=Math.max(0,k-(6-l));j<=Math.min(l,k);j++){
        const i=k-j;
        sum+=choose(6-l,i)*choose(l,j)*a**(6-l-i+j)*b**i*(-b)**(l-j);
      }
      out[k]+=c[l]*SQ(choose(6,l)/choose(6,k))*sum;
    }
    return out;
  }
  function filterNorm(k,c){return 1+8*k+(8+10*(1-c)*(1+3*c))*k*k+12*(1-c)*(1+c)**2*k**3;}
  function filterPoint(k,c){
    const Z=filterNorm(k,c),Zk=8+2*(8+10*(1-c)*(1+3*c))*k+36*(1-c)*(1+c)**2*k*k;
    const Zc=10*(2-6*c)*k*k+12*(1-2*c-3*c*c)*k**3;
    const x=k===0?0:2*(k*(1+k)*Zk*Zk+c*(1-c)*Zc*Zc/k)/(3*Z*Z);
    return [x,3*x-3+8/(Z*Z),2*x+8-8/(Z*Z)];
  }
  function fixedNormK(Z,c){
    let lo=0,hi=1;
    for(let i=0;i<1024&&filterNorm(hi,c)<Z;i++)hi*=2;
    if(!Number.isFinite(hi))return null;
    for(let i=0;i<70;i++){const mid=(lo+hi)/2;if(filterNorm(mid,c)<Z)lo=mid;else hi=mid;}
    return (lo+hi)/2;
  }
  function filteredPlane(point,tolerance=1e-7){
    if(!target(point,tolerance))return null;
    const [x,y,z]=point;
    if(y<5-tolerance||y>15+tolerance||Math.abs(z-(5*x+5-y))>tolerance)return null;
    const yy=clamp(y,5,15),t=SQ((yy-5)/10);
    if(x<(2/3)*t*(t+8)-tolerance||x>(yy+3)/3+tolerance)return null;
    let c=zero();
    if(near(point,[0,5,0],tolerance)){
      c[1]=c[5]=[1/SQ(2),0];
      return accepted(point,c,'Filtered-octahedral lower plane',{kind:'octahedral'},
        {latex:String.raw`\frac{D_1+D_5}{\sqrt2}`,parameterLatex:basisNote},tolerance);
    }
    const edge=[(yy+3)/3,yy,10+2*yy/3];
    if(near(point,edge,tolerance)){
      const p=1-SQ(Math.max((15-yy)/10,0));c[0]=[SQ(p),0];c[1]=[SQ(1-p),0];
      return accepted(point,c,'Product–W lower-plane edge',{p},
        {latex:String.raw`\sqrt p\,D_0+\sqrt{1-p}\,D_1`,parameterLatex:`p=1-\\sqrt{(15-S_2)/10},\\quad ${basisNote}`},tolerance);
    }
    const q=(y+3-3*x)/8;
    if(!(q>0&&q<1))return null;
    const Z=1/SQ(q);let lo=1/3,hi=1,k;
    for(let i=0;i<72;i++){
      const cc=(lo+hi)/2,kk=fixedNormK(Z,cc);if(kk===null)return null;
      if(filterPoint(kk,cc)[0]>x)lo=cc;else hi=cc;
    }
    const axis=(lo+hi)/2;k=fixedNormK(Z,axis);if(k===null)return null;
    const h=Math.asinh(SQ(k)),beta=Math.acos(SQ(axis));
    const rotated=rotateY([0,1/SQ(2),0,0,0,1/SQ(2),0],beta);
    // Removing the common exp(3h) prevents overflow near the limiting edge.
    c=rotated.map((a,k)=>[a*Math.exp(-h*k),0]);
    return accepted(point,c,'Filtered-octahedral lower plane',{k,c:axis,h,beta,Z,q},
      {latex:String.raw`\frac{e^{hJ_z}e^{-i\beta J_y}}{\sqrt Z}\,\frac{D_1+D_5}{\sqrt2}`,
        parameterLatex:`h\\approx ${number(h)},\\quad\\beta\\approx ${number(beta)},\\quad Z\\approx ${number(Z)},\\quad ${basisNote}`},tolerance);
  }

  function pairedOdd(point,tolerance){
    const [x,y]=point,f=SQ(1.5*Math.max(x,0)),w=5/8;
    const p=(1-w+f/2)/2,r=(1-w-f/2)/2,H=7.5*(y-5),h2=H-48*w*w;
    if(Math.min(p,r,h2)<-tolerance)return null;
    const A=SQ(Math.max(p,0)),C=SQ(Math.max(r,0)),den=2*A*C;
    let phase=0;
    if(den>1e-14){const cosine=(h2/(120*w)-(1-w))/den;if(Math.abs(cosine)>1+2*tolerance)return null;phase=Math.acos(clamp(cosine,-1,1))/2;}
    else if(Math.abs(h2-120*w*(1-w))>4*tolerance)return null;
    const c=zero();c[1]=[A,0];c[3]=[SQ(w)*Math.cos(phase),SQ(w)*Math.sin(phase)];c[5]=[C,0];
    return accepted(point,c,'Paired-purity face: odd C₂',{w,p1:p,p5:r,phase},
      {latex:String.raw`\sqrt{p_1}D_1+\sqrt{\frac58}e^{i\varphi}D_3+\sqrt{p_5}D_5`,
        parameterLatex:`p_1\\approx ${number(p)},\\quad p_5\\approx ${number(r)},\\quad\\varphi\\approx ${number(phase)},\\quad ${basisNote}`},tolerance);
  }
  function pairedEvenBridge(point,tolerance){
    const [x,,z]=point,u=1.5*Math.max(x,0),q=(8+2*x-z)/8;
    if(u>1/9+1e-12)return null;
    const rad=-240+2560*u+3840*q;
    if(rad<-2e-10)return null;
    const s=(68-SQ(Math.max(rad,0)))/128,f=SQ(u);
    if(s<1/4-1e-10||s>17/32+1e-10)return null;
    const populations=[(s+f/8)/2,(1-s+5*f/8)/2,(1-s-5*f/8)/2,(s-f/8)/2];
    if(Math.min(...populations)<-1e-10)return null;
    const [A,B,C,E]=populations.map(p=>SQ(Math.max(p,0))),R=4*A*A*E*E,W=4*B*B*C*C,den=2*SQ(R*W);
    let ell=0;
    if(den>1e-14){const cosine=(q-R-W)/den;if(Math.abs(cosine)>1+2e-9)return null;ell=Math.acos(clamp(cosine,-1,1));}
    else if(Math.abs(q-R-W)>1e-9)return null;
    const alpha=Math.atan2(C*E*Math.sin(ell),A*B+C*E*Math.cos(ell)),gamma=ell-alpha,c=zero();
    c[0]=[A*Math.cos(alpha),A*Math.sin(alpha)];c[2]=[B,0];c[4]=[C,0];c[6]=[E*Math.cos(gamma),E*Math.sin(gamma)];
    return accepted(point,c,'Paired-purity face: even C₂ bridge',{s,u,q,populations,ell,alpha,gamma},
      {latex:String.raw`\sqrt{p_0}e^{i\alpha}D_0+\sqrt{p_2}D_2+\sqrt{p_4}D_4+\sqrt{p_6}e^{i\gamma}D_6`,
        parameterLatex:String.raw`\begin{aligned}u&=3S_1/2,\quad f=\sqrt u,\quad q=(8+2S_1-S_3)/8,\\s&=\frac{68-\sqrt{-240+2560u+3840q}}{128},\\p_0&=(s+f/8)/2,\quad p_6=(s-f/8)/2,\\p_2&=(1-s+5f/8)/2,\quad p_4=(1-s-5f/8)/2,\\`+
          `\\alpha&\\approx ${number(alpha)},\\quad\\gamma\\approx ${number(gamma)},\\\\`+basisNote+String.raw`\end{aligned}`},tolerance);
  }
  function hodgeMatrix(n){
    const [a,b,c,d,e]=n;
    return [[0,0,0,0,SQ(10)*e/5,-d,c],
      [0,0,0,-2*SQ(5)*e/5,SQ(15)*d/5,0,-b],
      [0,0,2*SQ(6)*e/5,-SQ(2)*d/5,-3*c/5,SQ(15)*b/5,SQ(10)*a/5],
      [0,-2*SQ(5)*e/5,-SQ(2)*d/5,4*c/5,-SQ(2)*b/5,-2*SQ(5)*a/5,0],
      [SQ(10)*e/5,SQ(15)*d/5,-3*c/5,-SQ(2)*b/5,2*SQ(6)*a/5,0,0],
      [-d,0,SQ(15)*b/5,-2*SQ(5)*a/5,0,0,0],
      [c,-b,SQ(10)*a/5,0,0,0,0]];
  }
  function symmetricEigen(matrix){
    const n=matrix.length,A=matrix.map(r=>r.slice()),V=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>+(i===j)));
    for(let iteration=0;iteration<2000;iteration++){
      let p=0,q=1,max=0;
      for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)if(Math.abs(A[i][j])>max){max=Math.abs(A[i][j]);p=i;q=j;}
      if(max<2e-15)break;
      const tau=(A[q][q]-A[p][p])/(2*A[p][q]),t=(tau>=0?1:-1)/(Math.abs(tau)+SQ(1+tau*tau)),c=1/SQ(1+t*t),s=t*c,apq=A[p][q];
      A[p][p]-=t*apq;A[q][q]+=t*apq;A[p][q]=A[q][p]=0;
      for(let k=0;k<n;k++){
        if(k!==p&&k!==q){const a=A[k][p],b=A[k][q];A[k][p]=A[p][k]=c*a-s*b;A[k][q]=A[q][k]=s*a+c*b;}
        const a=V[k][p],b=V[k][q];V[k][p]=c*a-s*b;V[k][q]=s*a+c*b;
      }
    }
    return Array.from({length:n},(_,i)=>({value:A[i][i],vector:V.map(r=>r[i])}));
  }
  const hodgeLatex=String.raw`K(n)=\begin{pmatrix}0&0&0&0&\sqrt{10}n_4/5&-n_3&n_2\\0&0&0&-2\sqrt5n_4/5&\sqrt{15}n_3/5&0&-n_1\\0&0&2\sqrt6n_4/5&-\sqrt2n_3/5&-3n_2/5&\sqrt{15}n_1/5&\sqrt{10}n_0/5\\0&-2\sqrt5n_4/5&-\sqrt2n_3/5&4n_2/5&-\sqrt2n_1/5&-2\sqrt5n_0/5&0\\\sqrt{10}n_4/5&\sqrt{15}n_3/5&-3n_2/5&-\sqrt2n_1/5&2\sqrt6n_0/5&0&0\\-n_3&0&\sqrt{15}n_1/5&-2\sqrt5n_0/5&0&0&0\\n_2&-n_1&\sqrt{10}n_0/5&0&0&0&0\end{pmatrix}`;
  function pairedHodge(point,tolerance){
    const [x,,z]=point,q=(8+2*x-z)/8;
    if(x<2/27-1e-12||q<-1e-12)return null;
    const dd=(1+4*q-6*x)/5;if(dd<-1e-12)return null;
    const r=SQ(Math.max(q,0)),d=SQ(Math.max(dd,0));if(d>1+1e-12)return null;
    let r0,w,ell,branch;
    if(r>=d-1/4){
      const rad=1+4*r-4*d;if(rad<-1e-12)return null;
      r0=r;w=r-d;ell=(1-SQ(Math.max(rad,0)))/2;branch='first';
    }else{
      const rad=1+4*d-4*r;if(rad<-1e-12)return null;
      r0=-r;w=d-r;ell=(1+(3*d+2*r-1>=0?1:-1)*SQ(Math.max(rad,0)))/2;branch='second';
    }
    const denominator=r0+1-w;if(Math.abs(denominator)<1e-12)return null;
    const t=(ell-1)*(2*r0-w)/denominator,V=(r0+1-w)/2,U=1-V-(ell-t)**2/3-t*t;
    if(Math.min(U,V)<-2e-11)return null;
    const a=(ell-t)/SQ(6),n=[a,(SQ(Math.max(V,0))+SQ(Math.max(U,0)))/SQ(2),t,(SQ(Math.max(V,0))-SQ(Math.max(U,0)))/SQ(2),a];
    const eig=symmetricEigen(hodgeMatrix(n)).sort((a,b)=>Math.abs(a.value-1)-Math.abs(b.value-1));
    for(const pair of eig){
      if(Math.abs(pair.value-1)>2e-8)continue;
      const result=accepted(point,pair.vector.map(v=>[v,0]),'Paired-purity face: real Hodge bridge',{q,r,d,ell,t,U,V,n,branch,eigenvalue:pair.value},
        {latex:String.raw`|\psi\rangle=\sum_{k=0}^{6}c_kD_k,\qquad c\in\ker(K(n)-I_7),\quad\|c\|=1`,
          parameterLatex:`n\\approx\\left(${n.map(number).join(',\\,')}\\right),\\qquad ${hodgeLatex},\\qquad ${basisNote}`},tolerance);
      if(result)return result;
    }
    return null;
  }
  function pairedFace(point,tolerance=1e-7){
    if(!target(point,tolerance))return null;
    const [x,y,z]=point;if(x<-tolerance||x>3/8+tolerance||Math.abs(15+5*x-y-z)>tolerance)return null;
    return pairedOdd(point,tolerance)||pairedEvenBridge(point,tolerance)||pairedHodge(point,tolerance);
  }
  global.N6SurfaceConstructors=global.N6SurfaceConstructors||{};
  Object.assign(global.N6SurfaceConstructors,{filteredPlane,pairedFace});
})(typeof globalThis==='object'?globalThis:window);
