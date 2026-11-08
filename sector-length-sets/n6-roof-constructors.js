/* Explicit pure six-qubit q=0 roof inverses, ported from the Julia families.
 * Each returned state is normalized and independently checked through the
 * Dicke marginal-purity sector map. A failed inverse makes no exclusion claim.
 * No state-vector or sector-point interpolation is used for acceptance.
 */
(function (global) {
  'use strict';
  const api=global.N6SurfaceConstructors=global.N6SurfaceConstructors||{};
  const clamp=(x,a,b)=>Math.max(a,Math.min(b,x)),sqrt=x=>Math.sqrt(Math.max(0,x));
  const zero=()=>Array.from({length:7},()=>[0,0]);
  const add=(a,b)=>[a[0]+b[0],a[1]+b[1]],scale=(a,b)=>[a[0]*b,a[1]*b];
  const mul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]];
  const number=x=>Number(x.toPrecision(13));
  const discrepancy=(a,b)=>Math.max(...a.map((x,i)=>Math.abs(x-b[i])));
  function admissible(point,tolerance){return Array.isArray(point)&&point.length>=3&&point.slice(0,3).every(Number.isFinite)&&Number.isFinite(tolerance)&&tolerance>0&&Math.abs(point[2]-8-2*point[0])<=tolerance;}
  function checked(c,point,tolerance,family,parameters,latex,definitions,source){
    if(!global.SectorMath||!c||c.length!==7||c.some(z=>!z.every(Number.isFinite)))return null;
    const norm=c.reduce((s,z)=>s+z[0]*z[0]+z[1]*z[1],0);
    if(Math.abs(norm-1)>Math.max(tolerance,2e-12))return null;
    c=c.map(z=>scale(z,1/Math.sqrt(norm)));
    const sectors=global.SectorMath.sectorLengths(c);
    if(!sectors.every(Number.isFinite)||discrepancy(sectors.slice(0,3),point.slice(0,3))>tolerance)return null;
    return {N:6,status:'attained',coefficients:c,sectors,family,parameters,
      proofOrigin:source.split(';')[0],constructionSource:source,
      analyticalExpression:{latex,parameterLatex:definitions,
        analyticalStatus:'Explicit normalized pure-state family; recovered parameters checked by independent marginal purities'}};
  }
  function rank(point,tol){
    if(!admissible(point,tol))return null;
    let [x,y]=point;if(x < -tol||x>6+tol)return null;x=clamp(x,0,6);
    const f=sqrt(1.5*x),lower=Math.max(5+.4*(2-f)**2,3*x-3),upper=5+.4*(2+f)**2;
    if(y<lower-tol||y>upper+tol)return null;
    const target=clamp(y,lower,upper);
    const path=t=>{
      let d,e;
      if(f<=Math.sqrt(1.5)){
        const e0=sqrt(1-2*f*f/3),step=(1-e0)*Math.abs(t);
        e=e0+step;d=(-e+(t<0?-1:1)*sqrt(6*step*(2*e0+step)))/5;
      }else{
        const emin=f<=2?-1:1-2*sqrt(9-f*f)/Math.sqrt(5);
        e=emin+(1+t)*(1-emin)/2;d=(-e+sqrt(6*e*e+4*(f*f-1.5)))/5;
      }
      return [clamp(d,-1,1),clamp(e,-1,1)];
    };
    let a=-1,b=1;
    for(let i=0;i<90;i++){
      const t=(a+b)/2,[d,e]=path(t),ym=2.5*(3+2*d*d+2*d*e-d*d*e*e);
      if(Math.abs(ym-target)<2e-14){a=b=t;break;}
      if(ym<target)a=t;else b=t;
    }
    const [d,e]=path((a+b)/2),u=(1+d)/2,v=(1+e)/2,c=zero();
    c[0]=[sqrt(u*v),0];c[1]=[sqrt(u*(1-v)),0];c[5]=[sqrt((1-u)*v),0];c[6]=[sqrt((1-u)*(1-v)),0];
    return checked(c,point,tol,'Rank-one null-roof family',{u,v,d,e},
      '\\sqrt{uv}D_0+\\sqrt{u(1-v)}D_1+\\sqrt{(1-u)v}D_5+\\sqrt{(1-u)(1-v)}D_6',
      `u\\approx${number(u)},\\quad v\\approx${number(v)},\\qquad\\begin{aligned}d&=2u-1,\\quad e=2v-1,\\\\x&=\\frac{25d^2+10de-5e^2+6}{6},\\\\y&=\\frac52(3+2d^2+2de-d^2e^2).\\end{aligned}`,
      'Complete rank-one four-level image; Julia/N6_rankone_roof_attainment.jl');
  }
  function seven(point,tol){
    if(!admissible(point,tol))return null;
    const [x,y]=point;if(x<1/6-tol||x>8/3+tol)return null;
    const T=clamp((6*clamp(x,1/6,8/3)-1)/15,0,1);let delta=y-5;
    if(delta < -tol||delta>15*T*T/8+tol)return null;
    delta=clamp(delta,0,15*T*T/8);
    const denominator=T+sqrt(T*T-8*delta/15),t=clamp(denominator===0?0:(4*delta/15)/denominator,0,T);
    const lo=sqrt(1-T),hi=sqrt(T-t),r=(lo+2*Math.sqrt(2)*hi)/3,u=(2*Math.sqrt(2)*lo-hi)/3,v=sqrt(t);
    const a=sqrt((1+r)/2),b=[u/Math.sqrt(2*(1+r)),v/Math.sqrt(2*(1+r))],b2=mul(b,b);
    const q0=scale(add([a*a,0],b2),1/Math.sqrt(2)),q1=scale(add([a*a,0],scale(b2,-1)),1/Math.sqrt(2)),q2=scale(b,Math.sqrt(2)*a),c=zero();
    c[0]=add(scale(q0,Math.sqrt(2)/3),scale(q1,Math.sqrt(3)/4));
    c[2]=scale(q1,Math.sqrt(5)/4);c[3]=mul(q0,[0,Math.sqrt(5)/3]);
    c[4]=scale(q1,-Math.sqrt(5)/4);c[6]=add(scale(q0,Math.sqrt(2)/3),scale(q1,-Math.sqrt(3)/4));
    c[1]=mul(q2,[0,-1/Math.sqrt(2)]);c[5]=mul(q2,[0,1/Math.sqrt(2)]);
    return checked(c,point,tol,'Seven-level null-conic roof',{r,u,v,T,t},
      '\\frac{(a^2+b^2)|e_0\\rangle+(a^2-b^2)|e_1\\rangle+2ab|e_2\\rangle}{\\sqrt2}',
      `r\\approx${number(r)},\\quad u\\approx${number(u)},\\quad v\\approx${number(v)},\\qquad\\begin{aligned}a&=\\sqrt{(1+r)/2},\\quad b=\\frac{u+iv}{\\sqrt{2(1+r)}},\\\\|e_0\\rangle&=\\frac{\\sqrt2(D_0+D_6)+i\\sqrt5D_3}{3},\\\\|e_1\\rangle&=\\frac{\\sqrt3(D_0-D_6)+\\sqrt5(D_2-D_4)}4,\\\\|e_2\\rangle&=\\frac{i(D_5-D_1)}{\\sqrt2}.\\end{aligned}`,
      'Complete seven-level null-conic image; Julia/N6_sevenlevel_roof_attainment.jl');
  }
  function fixedDipoleCircle(point,tol){
    if(!admissible(point,tol)||!global.SectorMath)return null;
    const [x,y]=point;if(x<1/6-tol||x>8/3+tol)return null;
    const xx=clamp(x,1/6,8/3),f=sqrt(1.5*xx),ceiling=5+.4*(2-f)**2;
    if(y<5-tol||y>ceiling+tol)return null;
    const target=clamp(y,5,ceiling),T=clamp((6*xx-1)/15,0,1),ell=sqrt(1-T),h=sqrt(T);
    const r=(ell+2*Math.sqrt(2)*h)/3,u=(2*Math.sqrt(2)*ell-h)/3;
    const a=[Math.sqrt(6)*r/4,u,-Math.sqrt(10)*r/4],seed=[2/3,0,0,Math.sqrt(5)/3];
    const L=[[Math.sqrt(6)*u/2,-r/2,Math.sqrt(10)*u/2,-Math.sqrt(15)*r/2],
      [3*Math.sqrt(6)*r/4,2*u,-Math.sqrt(10)*r/4,0]];
    const dipole=L.map(row=>row.reduce((s,v,i)=>s+v*seed[i],0));
    const b0=seed.map((_,i)=>(L[0][i]*dipole[0]+L[1][i]*dipole[1])/4),w=seed.map((v,i)=>v-b0[i]);
    const v=[5*Math.sqrt(2)*(r*r-8*u*u)/96,
      5*r*(Math.sqrt(6)*r+4*Math.sqrt(3)*u)/48,
      Math.sqrt(30)*r*r/32+Math.sqrt(15)*r*u/6+Math.sqrt(30)*u*u/12,
      Math.sqrt(10)*(-r*r+8*u*u)/48];
    const state=theta=>{
      const b=b0.map((value,i)=>value+Math.cos(theta)*w[i]+Math.sin(theta)*v[i]),c=zero();
      // psi=(Eplus*a+Eminus*b)/sqrt2 in the real spin-flip bases.
      for(let k=0;k<3;k++){c[k]=[(a[k]+b[k])/2,0];c[6-k]=[((k%2?-1:1)*(a[k]-b[k]))/2,0];}
      c[3]=[b[3]/Math.sqrt(2),0];return c;
    };
    const secondSector=theta=>15*(4*global.SectorMath.reducedPurity(state(theta),2)-1-xx/3);
    // This is an analytical bracket, not an optimizer: y(0)=5, while
    // y(pi)=5+(4-f^2)(f^2+1)(17-3f^2)/50 >= ceiling throughout.
    let lo=0,hi=Math.PI,theta=0;
    if(Math.abs(target-5)>2e-14){
      if(secondSector(hi)<target-2e-12)return null;
      const desired=Math.min(2e-13,tol/16);
      for(let i=0;i<80;i++){
        theta=(lo+hi)/2;const yy=secondSector(theta);
        if(Math.abs(yy-target)<desired){lo=hi=theta;break;}
        if(yy<target)lo=theta;else hi=theta;
      }
      theta=(lo+hi)/2;
    }
    return checked(state(theta),point,tol,'Fixed-dipole null-roof circle',{x:xx,theta,r,u,a,b0,w,v},
      '\\frac1{\\sqrt2}\\left(\\sum_{k=0}^{2}a_k|E_k^+\\rangle+\\sum_{k=0}^{3}b_k(\\theta)|E_k^-\\rangle\\right)',
      `x=${number(xx)},\\quad\\theta\\approx${number(theta)},\\qquad\\begin{aligned}`+
      '|E_k^\\pm\\rangle&=\\frac{D_k\\pm(-1)^kD_{6-k}}{\\sqrt2}\\quad(0\\le k\\le2),\\\\|E_3^-\\rangle&=D_3,\\\\'+
      'T&=\\frac{6x-1}{15},\\quad\\ell=\\sqrt{1-T},\\quad h=\\sqrt T,\\\\r&=\\frac{\\ell+2\\sqrt2h}{3},\\quad u=\\frac{2\\sqrt2\\ell-h}{3},\\\\'+
      '\\mathbf a&=(\\sqrt6r/4,\\ u,\\ -\\sqrt{10}r/4),\\quad\\mathbf b_*=(2/3,0,0,\\sqrt5/3),\\\\'+
      'L&=\\begin{pmatrix}\\sqrt6u/2&-r/2&\\sqrt{10}u/2&-\\sqrt{15}r/2\\\\3\\sqrt6r/4&2u&-\\sqrt{10}r/4&0\\end{pmatrix},\\\\'+
      '\\mathbf b_0&=\\frac14L^TL\\mathbf b_*,\\quad\\mathbf w=\\mathbf b_*-\\mathbf b_0,\\\\'+
      '\\mathbf v&=\\left(\\frac{5\\sqrt2(r^2-8u^2)}{96},\\ \\frac{5r(\\sqrt6r+4\\sqrt3u)}{48},\\ \\frac{\\sqrt{30}r^2}{32}+\\frac{\\sqrt{15}ru}{6}+\\frac{\\sqrt{30}u^2}{12},\\ \\frac{\\sqrt{10}(-r^2+8u^2)}{48}\\right),\\\\'+
      '\\mathbf b(\\theta)&=\\mathbf b_0+\\cos\\theta\\,\\mathbf w+\\sin\\theta\\,\\mathbf v,\\qquad0\\le\\theta\\le\\pi,\\\\S_2(\\psi(x,\\theta))&=y.\\end{aligned}',
      'Complete fixed-dipole Stiefel-circle bridge; output/N6_threehour_border_search_20261020/roof/seven_edge_circle_proof.txt');
  }
  const binomial=[1,6,15,20,15,6,1];
  function coherent(r,phi){return binomial.map((b,k)=>{const amplitude=Math.sqrt(b)*Math.sqrt(r)**(6-k)*Math.sqrt(1-r)**k;return [amplitude*Math.cos(k*phi),amplitude*Math.sin(k*phi)];});}
  function reflected(point,tol,p,octa){
    const q=coherent(p.r,p.phi),c=q.map(z=>scale(z,-1));
    if(octa){c[1]=q[5];c[5]=q[1];}else{c[0]=q[6];c[6]=q[0];}
    const name=octa?'o':'g',definition=octa?'\\frac{D_1+D_5}{\\sqrt2}':'\\frac{D_0+D_6}{\\sqrt2}';
    return checked(c,point,tol,octa?'Octa–product roof':'GHZ–product roof',p,
      `(2|${name}\\rangle\\langle ${name}|-I)|p(r,\\phi)\\rangle`,
      `r\\approx${number(p.r)},\\quad\\phi\\approx${number(p.phi)},\\qquad\\begin{aligned}|${name}\\rangle&=${definition},\\\\|p(r,\\phi)\\rangle&=(\\sqrt r|0\\rangle+e^{i\\phi}\\sqrt{1-r}|1\\rangle)^{\\otimes6}.\\end{aligned}`,
      'Complete Householder-reflected coherent-state image; Julia/N6_'+(octa?'octa':'GHZ')+'_product_attainment.jl');
  }
  function gp(point,tol){
    if(!admissible(point,tol))return null;
    const [x,y,z]=point;if(discrepancy(point,[6,15,20])<=tol)return reflected(point,tol,{r:1,phi:0,t:0,C:1,L:0},false);
    const X=6-x,Y=15-y,L=(3*X+Y)/4;if(!(L>0))return null;
    const map=(t,C)=>{const l=120*t*((1-3*t)**2-2*(1-C)*t**3),xx=6-(1-t)*l;return [xx,15-(1+3*t)*l,8+2*xx];};
    if(L<1e-3){
      let t=L/120;for(let i=0;i<4;i++)t-=(120*t*(1-6*t+9*t*t-2*t**3)-L)/(120*(1-12*t+27*t*t-8*t**3));
      if(t>0&&t<.25&&discrepancy(map(t,0),point)<=tol)return reflected(point,tol,{r:2*t/(1+sqrt(1-4*t)),phi:Math.PI/12,t,C:0,L},false);
    }
    let t=(Y-X)/(4*L);const slack=tol/Math.max(L,Number.EPSILON);
    if(t < -slack||t>.25+slack)return null;t=clamp(t,0,.25);
    const min=120*t*(1-t)**2*(1-4*t),max=120*t*(1-3*t)**2;
    if(L<min-tol||L>max+tol||t===0)return null;
    const width=480*t**4,C=width===0?0:clamp(1-2*(max-L)/width,-1,1);
    if(discrepancy(map(t,C),point)>tol)return null;
    return reflected(point,tol,{r:2*t/(1+sqrt(1-4*t)),phi:Math.acos(C)/6,t,C,L},false);
  }
  function opInvariants(E,F,point,tol){
    const slack=Math.max(512*Number.EPSILON,8*tol);
    if(E < -slack||E>1/3+slack||F < -slack||F>1/27+slack)return null;
    E=clamp(E,0,1/3);F=clamp(F,0,1/27);
    if(E*E-4*E**3-4*F-27*F*F+18*E*F < -slack)return null;
    const h=Math.max(1-3*E,0);let axes;
    if(F===0){if(E>.25+slack)return null;const small=2*E/(1+sqrt(1-4*E));axes=[0,small,1-small];}
    else if(h<=256*Number.EPSILON)axes=[1/3,1/3,1/3];
    else{
      const angle=Math.acos(clamp((2-9*E+27*F)/(2*h*Math.sqrt(h)),-1,1))/3,amplitude=2*Math.sqrt(h)/3;
      axes=[0,1,2].map(k=>1/3+amplitude*Math.cos(angle-2*k*Math.PI/3));
      if(Math.min(...axes)<-slack)return null;axes=axes.map(x=>Math.max(x,0));const total=axes.reduce((a,b)=>a+b,0);axes=axes.map(x=>x/total);
    }
    const largest=axes.indexOf(Math.max(...axes)),others=axes.filter((_,i)=>i!==largest),[A,B]=others;
    const t=(A+B)/4,C=t===0?1:clamp(8*A*B/(A+B)**2-1,-1,1);
    const ee=4*t+(2*C-14)*t*t,ff=2*t*t*(1-4*t)*(1+C),D=ee-ff,xx=6-30*D*(1-ee),actual=[xx,15-30*D*(1+3*ff),8+2*xx];
    if(discrepancy(actual,point)>tol)return null;
    return reflected(point,tol,{r:(1+sqrt(1-4*t))/2,phi:Math.acos(C)/4,t,C,E,F,axisSquares:[A,B,axes[largest]]},true);
  }
  function op(point,tol){
    if(!admissible(point,tol))return null;
    const [x,y]=point;if(discrepancy(point,[6,15,20])<=tol)return reflected(point,tol,{r:1,phi:0,t:0,C:1,E:0,F:0},true);
    const X=(6-x)/30,Y=(15-y)/30,U=3*X+Y;if(!(U>0)||4-3*U<0)return null;
    const D=U/(2+sqrt(4-3*U));if(!(D>0)||D>8/27+tol)return null;
    if(D<1e-3){const candidate=opInvariants(D,0,point,tol);if(candidate)return candidate;}
    const F=(Y/D-1)/3;return opInvariants(D+F,F,point,tol);
  }
  function ww(point,tol){
    if(!admissible(point,tol))return null;
    const dx=8/3-point[0],dy=point[1]-5;let A=18*(163*dy-96*dx)/1625,B=18*(21*dx-32*dy)/325;
    if(Math.min(A,B)<-tol)return null;A=Math.max(A,0);B=Math.max(B,0);
    const k=1+2*Math.sqrt(A)/9,delta=k*k-4*(A+B);if(delta < -tol)return null;
    const N=2/(k+sqrt(delta)),a=N*Math.sqrt(A),b=N*Math.sqrt(B),norm=Math.sqrt(1+a*a+b*b-2*a/9);
    const vector=[-Math.sqrt(30)/216,-1/9,-5*Math.sqrt(2)/24,-5*Math.sqrt(30)/54,-25*Math.sqrt(2)/72,0,25*Math.sqrt(30)/216];
    const c=vector.map((v,i)=>[((i===1?1:0)+a*v)/norm,b*v/norm]);
    return checked(c,point,tol,'Rotated W–W roof',{a,b,A,B,N},
      '\\frac{|u\\rangle+(a+ib)|v\\rangle}{\\sqrt{1+a^2+b^2-2a/9}}',
      `a\\approx${number(a)},\\quad b\\approx${number(b)},\\qquad\\begin{aligned}|u\\rangle&=D_1,\\\\|v\\rangle&=-\\frac{\\sqrt{30}}{216}D_0-\\frac19D_1-\\frac{5\\sqrt2}{24}D_2-\\frac{5\\sqrt{30}}{54}D_3-\\frac{25\\sqrt2}{72}D_4+\\frac{25\\sqrt{30}}{216}D_6.\\end{aligned}`,
      'Complete two-rotated-W image; Julia/N6_rotated_W_attainment.jl');
  }
  function newton2(point,seed,map,{iterations=45,determinantTolerance=1e-15,maxStep=[.15,.20],inDomain=()=>true,tolerance=1e-7}={}){
    let [a,b]=seed;const desired=Math.min(2e-13,tolerance/16);
    for(let i=0;i<iterations;i++){
      const [x,y,xa,xb,ya,yb]=map(a,b),rx=x-point[0],ry=y-point[1],error=Math.max(Math.abs(rx),Math.abs(ry));
      if(error<=desired)return [a,b];
      const determinant=xa*yb-xb*ya;if(Math.abs(determinant)<=determinantTolerance)return null;
      let da=(rx*yb-ry*xb)/determinant,db=(xa*ry-ya*rx)/determinant;
      const scale=Math.max(1,Math.abs(da)/maxStep[0],Math.abs(db)/maxStep[1]);da/=scale;db/=scale;
      let accepted=false,step=1;
      for(let j=0;j<16;j++,step/=2){
        const aa=a-step*da,bb=b-step*db;
        if(inDomain(aa,bb)){
          const [xx,yy]=map(aa,bb);
          if(Math.max(Math.abs(xx-point[0]),Math.abs(yy-point[1]))<error){a=aa;b=bb;accepted=true;break;}
        }
      }
      if(!accepted)return null;
    }
    return null;
  }
  function smallFourMap(t,k){
    const f=t*(2+3*k)-2,h=1+9*k*k;
    return [(2/3)*f*f+t*(1-t)*(1-k)**2,
      5+t*t*h*h/10+t*(1-t)*(13-10*k-3*k*k),
      (4/3)*f*(2+3*k)+(1-2*t)*(1-k)**2,
      4*t*f-2*t*(1-t)*(1-k),
      t*h*h/5+(1-2*t)*(13-10*k-3*k*k),
      (18/5)*t*t*k*h+t*(1-t)*(-10-6*k)];
  }
  const smallFourSeeds=[[5/8,1/3],[313/500,1/3],[.5,1]];
  for(let i=0;i<36;i++)for(let j=0;j<37;j++){
    const t=.30+.70*i/35,k=-.20+1.20*j/36,[x,y]=smallFourMap(t,k);
    if(x>=0&&x<=.19&&y>=5&&y<=7.6)smallFourSeeds.push([t,k]);
  }
  function smallfour(point,tol){
    if(!admissible(point,tol)||point[0]<-tol||point[0]>1/6+tol||point[1]<5-tol||point[1]>7.5+tol)return null;
    const seeds=smallFourSeeds.map(seed=>{const [x,y]=smallFourMap(...seed);return {seed,d:(x-point[0])**2+(y-point[1])**2/4};}).sort((a,b)=>a.d-b.d).slice(0,10);
    for(const {seed}of seeds){
      const fit=newton2(point,seed,smallFourMap,{tolerance:tol,inDomain:(t,k)=>t>=0&&t<=1&&k>=-1&&k<=1});if(!fit)continue;
      const [t,k]=fit,c=zero();c[0]=[Math.sqrt(t)*(1+k)/2,0];c[3]=[sqrt(t*(1-k*k)/2),0];c[5]=[-sqrt(1-t),0];c[6]=[Math.sqrt(t)*(1-k)/2,0];
      const result=checked(c,point,tol,'Four-level null roof',{t,k},
        '\\frac{\\sqrt t(1+k)}2D_0+\\sqrt{\\frac{t(1-k^2)}2}D_3-\\sqrt{1-t}D_5+\\frac{\\sqrt t(1-k)}2D_6',
        `t\\approx${number(t)},\\quad k\\approx${number(k)},\\qquad0\\le t\\le1,\\quad-1\\le k\\le1.`,
        'Target-specific four-level null reconstruction; Julia/N6_smallx_fourlevel_attainment.jl');
      if(result)return result;
    }
    return null;
  }
  function fiveMap(p,a,b){
    const ca=Math.cos(a),eta=Math.cos(2*a),sa2=Math.sin(2*a),cb=Math.cos(b),sb=Math.sin(b),nu=Math.cos(2*b),sb2=Math.sin(2*b);
    const f=2*eta+nu,fb=-2*sb2,q=1+4*eta*nu,qb=-8*eta*sb2,K=Math.sqrt(60)*eta+6*sa2,t=sb2*K/2,tb=nu*K;
    const h=ca*(Math.sqrt(6)*cb-Math.sqrt(10)*sb),hb=ca*(-Math.sqrt(6)*sb-Math.sqrt(10)*cb);
    const g=ca*(5*Math.sqrt(6)*cb-3*Math.sqrt(10)*sb)/2,gb=ca*(-5*Math.sqrt(6)*sb-3*Math.sqrt(10)*cb)/2,w=1-p,F=2*p+w*f;
    return [(2/3)*(F*F+p*w*h*h),
      5+(2/5)*w*w*q*q+(8/15)*p*w*g*g+(2/15)*w*w*t*t,
      (2/3)*(2*F*(2-f)+(1-2*p)*h*h),
      (2/3)*(2*F*w*fb+2*p*w*h*hb),
      -(4/5)*w*q*q+(8/15)*(1-2*p)*g*g-(4/15)*w*t*t,
      (4/5)*w*w*q*qb+(16/15)*p*w*g*gb+(4/15)*w*w*t*tb];
  }
  // Initial guesses only, copied from the existing five-level Julia inverse.
  // Every proposed fit is solved at the clicked target and physically checked.
  const fiveSeeds=[[0.5,-0.5235987755982988,1.0471975511965976],[0.1549613889652307,-0.48129576455512074,0.9742325183313576],[0.1764813741432611,-0.4793134629882564,0.9767534154786334],[0.19603706539557647,-0.4794184862558651,0.9701550504713263],[0.18197753097131225,-0.4502866966695081,-2.1496679149409665],[0.23928805243386544,-0.4791869842631192,0.9646688157620729],[0.35258562090796103,-0.5458436749148485,-2.167817953850956],[0.26220904605620726,-0.478569984277361,0.9659006832927773],[0.7676197900259607,1.3992997570379295,0.5357433311541253],[0.15038754169007668,-0.413969873061967,-2.1967617594477113],[0.7629193706153563,1.2582462116200948,0.5727876292702017],[0.31046885053622664,-0.48159705627899474,0.9588343984084612],[0.8086460285552097,-1.5563538719455758,-2.4484276527627475],[0.7749247833714799,-1.5895051044134718,-0.358863468243698],[0.21763980008293538,-0.4313322652388298,0.9245809204614726],[0.812920723329814,1.1758976104002876,-2.1370459237783237],[0.7813367782550438,-1.5848875565709428,-0.3276443874129322],[0.36141532150954175,2.6566501821121222,-2.1849044733837144],[0.7879415907870297,1.2592588453256957,0.5529623826485806],[0.8295851961125705,-1.478404854927221,-0.7963307052078482],[0.8644121801867028,1.5661804019692367,-1.427791242622842],[0.3789580280839747,-0.4797589256771835,0.9649681159507522],[0.8104091636482997,1.5290449166950453,-0.4717150257491324],[0.7371630305223369,0.9331806827730024,0.7230741129914214],[0.34583662450774477,-0.44367289879959254,0.9749562708775195],[0.790557463330256,1.1183466185813953,0.6388433808610369],[0.8001561691748377,1.5559688831389111,0.21233103535525552],[0.828105725622813,1.1472117723202226,0.8554530221476583],[0.8613099518686945,1.5148978457087052,-0.929052968382824],[0.8748362123333974,-1.50855706669226,-1.2083548961592838],[0.8455368087991247,-1.7992162660147402,0.7132743136622919],[0.8728092807587178,1.574825484810728,-1.009486499126392],[0.8302782265703317,-2.1129135555210197,-1.9937874589008944],[0.5026404371751562,-0.5066083206562508,-2.198770141055174],[0.38859976978473115,-0.4195054851488881,1.013340422876938],[0.8304891299250832,1.0311422170630773,0.8458470841194954],[0.48196444866749655,-0.46616133572565494,1.0197389076187888],[0.3621213019488308,-0.39734929550187853,0.9298262124077862],[0.8594278635538802,-1.4811879633601481,-0.5419221710288074],[0.5921099478365256,-0.5285082839105741,0.9317080170752805],[0.6338859863079522,-0.5830556167408169,0.856195572545791],[0.8705216510602021,-1.8762211373605573,0.6592022583191554],[0.7864384837778735,0.775637783347343,0.8844565866261536],[0.8455985144204029,0.9267336838788741,0.8920988072627606],[0.8994992720912434,-1.5181054665424325,0.7383098706871553],[0.5967620591747564,-0.47444680995307903,0.8833761510494351],[0.7502565497306073,2.541549146734173,-2.208159832705066],[0.927476391537395,1.558628757418839,-1.0211680145608852],[0.8194145561353452,-0.6730265137192264,0.9247435385886954],[0.6086627239346084,-0.4116241792010213,-2.26192106421967],[0.9400087619159876,-1.3154065802510415,-2.0390530866669123]];
  function five(point,tol){
    if(!admissible(point,tol)||point[0]<.15-tol||point[0]>2.7+tol||point[1]<5-tol||point[1]>5.25+tol)return null;
    const seeds=fiveSeeds.map(seed=>{const [x,y]=fiveMap(...seed);return {seed,d:(x-point[0])**2+64*(y-point[1])**2};}).sort((a,b)=>a.d-b.d).slice(0,8);
    for(const {seed:[pp,a,bb]}of seeds){
      const map=(p,b)=>fiveMap(p,a,b),fit=newton2(point,[pp,bb],map,{iterations:35,determinantTolerance:1e-13,maxStep:[.20,.30],inDomain:p=>p>=0&&p<=1,tolerance:tol});if(!fit)continue;
      const [p,b]=fit,ca=Math.cos(a),sa=Math.sin(a),cb=Math.cos(b),sb=Math.sin(b),scale=Math.sqrt(1-p),c=zero();
      c[0]=[scale*ca*cb,0];c[1]=[0,Math.sqrt(p)];c[2]=[scale*ca*sb,0];c[4]=[scale*sa*cb,0];c[6]=[-scale*sa*sb,0];
      const result=checked(c,point,tol,'Five-level null roof',{p,a,b},
        '\\sqrt{1-p}\\left(\\cos a\\,\\cos b\\,D_0+\\cos a\\,\\sin b\\,D_2+\\sin a\\,\\cos b\\,D_4-\\sin a\\,\\sin b\\,D_6\\right)+i\\sqrt p\\,D_1',
        `p\\approx${number(p)},\\quad a\\approx${number(a)},\\quad b\\approx${number(b)},\\qquad0\\le p\\le1.`,
        'Target-specific five-level null reconstruction; Julia/N6_fivelevel_roof_attainment.jl');
      if(result)return result;
    }
    return null;
  }
  const families={rank,seven,fixedDipoleCircle,gp,op,ww,smallfour,five};
  api.roof=function(point,tolerance=1e-7){
    if(!admissible(point,tolerance))return null;
    for(const construct of Object.values(families)){const result=construct(point,tolerance);if(result)return result;}
    return null;
  };
  api.roofFamilies=Object.freeze(families);
})(typeof window==='object'?window:globalThis);
