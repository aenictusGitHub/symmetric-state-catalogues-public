(function (root) {
  'use strict';
  const choose2 = n => n*(n-1)/2;

  function fromSectors(n, [s1, s2]) {
    const p1 = (1+s1/n)/2;
    if (n <= 3) return [p1];
    return [p1, (1+2*s1/n+s2/choose2(n))/4];
  }

  function toSectors(n, [p1, p2]) {
    const s1 = n*(2*p1-1);
    if (n <= 3) return [s1];
    return [s1, choose2(n)*(4*p2-1-2*s1/n)];
  }

  function exactPurityLabels(n, sectors) {
    const parse = value => {
      if (!/^-?\d+(?:\/[1-9]\d*)?$/.test(value ?? '')) return null;
      const [a,b='1'] = value.split('/');
      return [Number(a),Number(b)];
    };
    const fraction = (a,b) => {
      let x=Math.abs(a),y=b;
      while(y) [x,y]=[y,x%y];
      return b===x ? String(a/x) : '\\frac{'+a/x+'}{'+b/x+'}';
    };
    const x=parse(sectors?.[0]), y=parse(sectors?.[1]);
    const p1=x ? fraction(n*x[1]+x[0],2*n*x[1]) : null;
    const c=choose2(n);
    const p2=n===3?p1:x&&y?fraction(n*c*x[1]*y[1]+2*c*x[0]*y[1]+n*x[1]*y[0],
                                      4*n*c*x[1]*y[1]):null;
    return Array.from({length:n},(_,k)=>{
      const t=Math.min(k+1,n-k-1);
      return t===0?'1':t===1?p1:p2;
    });
  }

  const lineAxes = {ranges:[[.5,1],[0,1]], ticks:[[.5,.6,.7,.8,.9,1],[]]};
  const axes = {
    2: lineAxes,
    3: lineAxes,
    4: {ranges:[[.5,1],[1/3,1]], ticks:[[.5,.6,.7,.8,.9,1],[1/3,.5,.6,.7,.8,.9,1]]},
    5: {ranges:[[.5,1],[1/3,1]], ticks:[[.5,.6,.7,.8,.9,1],[1/3,.5,.6,.7,.8,.9,1]]}
  };
  const api = {fromSectors,toSectors,exactPurityLabels,axes};
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.PurityCoordinates = api;
})(globalThis);
