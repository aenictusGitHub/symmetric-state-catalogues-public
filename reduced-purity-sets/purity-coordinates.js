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

  const lineAxes = {ranges:[[.5,1],[0,1]], ticks:[[.5,.6,.7,.8,.9,1],[]]};
  const axes = {
    2: lineAxes,
    3: lineAxes,
    4: {ranges:[[.5,1],[1/3,1]], ticks:[[.5,.6,.7,.8,.9,1],[1/3,.5,.6,.7,.8,.9,1]]},
    5: {ranges:[[.5,1],[1/3,1]], ticks:[[.5,.6,.7,.8,.9,1],[1/3,.5,.6,.7,.8,.9,1]]}
  };
  const api = {fromSectors,toSectors,axes};
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.PurityCoordinates = api;
})(globalThis);
