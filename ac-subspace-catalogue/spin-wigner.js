/* Unit-integral spin Wigner function for the catalogue's normalized rho=P_V/d.
 *
 * The all-positive Stratonovich-Weyl kernel is
 *   Delta(n)=sqrt(4*pi/(N+1)) sum_lm T_lm Y_lm(n)* .
 * See Sanchez-Soto et al., Phase Space Insights: Wigner Functions for Qubits and
 * Beyond (2025), https://doi.org/10.3390/app15095155, Sec. 2.
 * We multiply Tr(rho Delta) by (N+1)/(4*pi), so integral W dOmega=1.
 * In the existing OperatorMath convention rho_lm=Tr(rho T_lm^dagger),
 * ordered m=l,...,-l, this is
 *   W(n)=sqrt((N+1)/(4*pi)) sum_lm rho_lm Y_lm(n).
 * The absence of conjugation here is intentional; the kernel uses conjugate
 * harmonics, whereas its trace contains Tr(rho T_lm). Hermiticity converts
 * that trace to the coefficient formula above. N=1 gives
 *   W(n)=(1+sqrt(3) r dot n)/(4*pi).
 *
 * Orthonormal Condon-Shortley harmonics are evaluated with normalized
 * associated-Legendre recurrences (no large factorials). A regular grid
 * reuses its polar Fourier amplitudes and azimuthal phases.
 */
(function(root){
  'use strict';
  const PI4=4*Math.PI, BASE=1/Math.sqrt(PI4);
  const reference='https://doi.org/10.3390/app15095155';
  const complex=z=>Array.isArray(z)?[Number(z[0]),Number(z[1]||0)]:[Number(z),0];
  function prepare(model){
    const n=Number(model?.n);
    if(!Number.isInteger(n)||n<0||!Array.isArray(model?.multipoles))throw new Error('A finite spin-multipole model is required.');
    const byRank=new Map(model.multipoles.map(b=>[Number(b.rank),b]));
    const re=Array.from({length:n+1},()=>new Float64Array(n+1)),im=Array.from({length:n+1},()=>new Float64Array(n+1));
    let hermiticityResidual=0,coefficientNorm=0;
    for(let l=0;l<=n;l++){
      const block=byRank.get(l);
      if(!block||!Array.isArray(block.coefficients)||block.coefficients.length!==2*l+1)throw new Error('Complete multipole coefficients are required through rank N.');
      for(let m=0;m<=l;m++){
        const positive=complex(block.coefficients[l-m]),negative=complex(block.coefficients[l+m]),phase=m%2?-1:1;
        if(!positive.every(Number.isFinite)||!negative.every(Number.isFinite))throw new Error('Multipole coefficients must be finite.');
        hermiticityResidual=Math.max(hermiticityResidual,Math.hypot(negative[0]-phase*positive[0],negative[1]+phase*positive[1]));
        coefficientNorm=Math.max(coefficientNorm,Math.hypot(...positive));re[m][l]=positive[0];im[m][l]=positive[1];
      }
    }
    if(hermiticityResidual>1e-8*Math.max(1,coefficientNorm))throw new Error('A real spin Wigner function requires Hermitian multipole coefficients.');
    const factor=Math.sqrt((n+1)/PI4),aa=Array.from({length:n+1},()=>new Float64Array(n+1)),bb=Array.from({length:n+1},()=>new Float64Array(n+1));
    for(let m=0;m<=n;m++)for(let l=m+2;l<=n;l++){
      aa[m][l]=Math.sqrt((4*l*l-1)/(l*l-m*m));
      bb[m][l]=Math.sqrt((2*l+1)*((l-1)*(l-1)-m*m)/((2*l-3)*(l*l-m*m)));
    }
    function fourier(z){
      const c=new Float64Array(n+1),s=new Float64Array(n+1),sin=Math.sqrt(Math.max(0,1-z*z));let diagonal=BASE;
      for(let m=0;m<=n;m++){
        if(m)diagonal*=-Math.sqrt((2*m+1)/(2*m))*sin;
        let previous=diagonal;c[m]+=previous*re[m][m];s[m]+=previous*im[m][m];
        if(m===n)continue;
        let current=z*Math.sqrt(2*m+3)*previous;c[m]+=current*re[m][m+1];s[m]+=current*im[m][m+1];
        for(let l=m+2;l<=n;l++){
          const next=aa[m][l]*z*current-bb[m][l]*previous;c[m]+=next*re[m][l];s[m]+=next*im[m][l];previous=current;current=next;
        }
      }
      return {c,s};
    }
    function azimuth(amplitudes,phi){
      const {c,s}=amplitudes,cp=Math.cos(phi),sp=Math.sin(phi);let cos=1,sin=0,sum=c[0];
      for(let m=1;m<=n;m++){const next=cos*cp-sin*sp;sin=sin*cp+cos*sp;cos=next;sum+=2*(c[m]*cos-s[m]*sin);}
      return factor*sum;
    }
    function at(theta,phi){
      theta=Number(theta);phi=Number(phi);if(!Number.isFinite(theta)||!Number.isFinite(phi))throw new Error('Finite spherical angles are required.');
      return azimuth(fourier(Math.max(-1,Math.min(1,Math.cos(theta)))),phi);
    }
    function value(direction){
      const v=Array.isArray(direction)||ArrayBuffer.isView(direction)?direction:[direction?.x,direction?.y,direction?.z],x=Number(v[0]),y=Number(v[1]),z=Number(v[2]),length=Math.hypot(x,y,z);
      if(!(length>0)||![x,y,z].every(Number.isFinite))throw new Error('A nonzero finite spatial direction is required.');
      return azimuth(fourier(Math.max(-1,Math.min(1,z/length))),Math.atan2(y,x));
    }
    function grid(rows=128,cols=256,options={}){
      rows=Number(rows);cols=Number(cols);
      if(!Number.isInteger(rows)||!Number.isInteger(cols)||rows<1||cols<1)throw new Error('Positive integer grid dimensions are required.');
      const values=new Float64Array((rows+1)*(cols+1)),cosines=Array.from({length:cols},()=>new Float64Array(n+1)),sines=Array.from({length:cols},()=>new Float64Array(n+1));
      for(let j=0;j<cols;j++){
        const phi=2*Math.PI*j/cols,cp=Math.cos(phi),sp=Math.sin(phi);let cos=1,sin=0;cosines[j][0]=1;
        for(let m=1;m<=n;m++){const next=cos*cp-sin*sp;sin=sin*cp+cos*sp;cos=next;cosines[j][m]=cos;sines[j][m]=sin;}
      }
      let min=Infinity,max=-Infinity;
      for(let i=0;i<=rows;i++){
        const {c,s}=fourier(Math.cos(Math.PI*i/rows)),offset=i*(cols+1);
        for(let j=0;j<cols;j++){
          let sum=c[0];for(let m=1;m<=n;m++)sum+=2*(c[m]*cosines[j][m]-s[m]*sines[j][m]);
          const w=factor*sum;values[offset+j]=w;min=Math.min(min,w);max=Math.max(max,w);
        }
        values[offset+cols]=values[offset];
      }
      if(options.geometry===false)return {rows,cols,values,min,max,qMin:min,qMax:max};
      const vertices=[],faces=[],directions=Array.from({length:cols},(_,j)=>[Math.cos(2*Math.PI*j/cols),Math.sin(2*Math.PI*j/cols)]);
      for(let i=0;i<=rows;i++){
        const theta=Math.PI*i/rows,z=Math.cos(theta),radius=i===0||i===rows?0:Math.sin(theta);
        for(let j=0;j<=cols;j++){const [x,y]=directions[j===cols?0:j];vertices.push({n:[radius*x,radius*y,z],q:values[i*(cols+1)+j]});}
      }
      for(let i=0;i<rows;i++)for(let j=0;j<cols;j++){const a=i*(cols+1)+j,b=a+cols+1;faces.push([a,a+1,b+1,b]);}
      return {rows,cols,values,vertices,faces,min,max,qMin:min,qMax:max};
    }
    return {n,value,at,grid,normalization:'integral W dOmega = 1',formula:'W(n)=sqrt((N+1)/(4*pi))*sum_lm rho_lm Y_lm(n)',isotropicValue:1/PI4,hermiticityResidual,reference};
  }
  const api={prepare,reference,version:'1.0.0'};root.SpinWigner=api;if(typeof module==='object'&&module.exports)module.exports=api;
})(typeof window==='object'?window:globalThis);
