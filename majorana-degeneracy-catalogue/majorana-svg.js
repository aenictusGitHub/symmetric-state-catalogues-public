/* Vector-only SVG recorder for the Majorana renderer's Canvas2D drawing subset.
 * Keeping a single drawing routine preserves view orientation and visible layers.
 * All SVG content is local; labels/metadata are XML-escaped. */
(() => {
  'use strict';
  const xml=value=>String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&apos;'}[c]));
  const number=value=>{if(!Number.isFinite(value))throw new Error('Non-finite SVG coordinate');return String(Number(value.toFixed(6)));};
  // Inkscape requires alpha as an SVG opacity attribute, not rgba() paint.
  const paint=(value,attribute,opacityAttribute=attribute+'-opacity')=>{
    if(value?.id)return `${attribute}="url(#${value.id})"`;
    const rgba=/^rgba\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*\)$/i.exec(String(value));
    const color=rgba?`rgb(${rgba.slice(1,4).map(v=>Math.max(0,Math.min(255,Math.round(Number(v))))).join(',')})`:value==='transparent'?'#000':value;
    const alpha=rgba?Math.max(0,Math.min(1,Number(rgba[4]))):value==='transparent'?0:1;
    return `${attribute}="${xml(color)}" ${opacityAttribute}="${number(alpha)}"`;
  };
  class Context {
    constructor(){
      this.elements=[];this.defs=[];this.stack=[];this.path=[];
      this.state={fillStyle:'#000',strokeStyle:'#000',globalAlpha:1,lineWidth:1,lineJoin:'miter',font:'12px sans-serif',textAlign:'left',textBaseline:'alphabetic',shadowColor:'transparent',shadowBlur:0,shadowOffsetX:0,shadowOffsetY:0,lineDashOffset:0,dash:[],transform:[1,0,0,1,0,0]};
      for(const key of Object.keys(this.state))Object.defineProperty(this,key,{get:()=>this.state[key],set:value=>{this.state[key]=value;}});
    }
    save(){this.stack.push({...this.state});}
    restore(){if(this.stack.length)this.state=this.stack.pop();}
    setTransform(...values){this.state.transform=values;}
    translate(x,y){const [a,b,c,d,e,f]=this.state.transform;this.state.transform=[a,b,c,d,e+a*x+c*y,f+b*x+d*y];}
    scale(x,y){const [a,b,c,d,e,f]=this.state.transform;this.state.transform=[a*x,b*x,c*y,d*y,e,f];}
    setLineDash(value){this.state.dash=[...value];}
    clearRect(){/* Each export starts with an empty recording. */}
    beginPath(){this.path=[];}
    moveTo(x,y){this.path.push(`M${number(x)} ${number(y)}`);}
    lineTo(x,y){this.path.push(`L${number(x)} ${number(y)}`);}
    closePath(){this.path.push('Z');}
    arc(x,y,r,start,end){
      // All circles used by this renderer are complete and start a new subpath.
      if(Math.abs(end-start)<2*Math.PI-1e-9)throw new Error('SVG recorder expects complete circles');
      const n=number;this.path.push(`M${n(x+r)} ${n(y)} A${n(r)} ${n(r)} 0 1 0 ${n(x-r)} ${n(y)} A${n(r)} ${n(r)} 0 1 0 ${n(x+r)} ${n(y)} Z`);
    }
    createRadialGradient(x0,y0,r0,x1,y1,r1){
      const id=`gradient-${this.defs.length}`,stops=[],start=r0/r1;
      // SVG 1.1 has no focal-radius attribute. Extrapolate the zero-radius
      // focus and remap stops to reproduce the canvas's two nested circles.
      const fx=(x0-start*x1)/(1-start),fy=(y0-start*y1)/(1-start);
      this.defs.push(()=>`<radialGradient id="${id}" gradientUnits="userSpaceOnUse" cx="${number(x1)}" cy="${number(y1)}" r="${number(r1)}" fx="${number(fx)}" fy="${number(fy)}">${stops.join('')}</radialGradient>`);
      return {id,addColorStop(offset,color){stops.push(`<stop offset="${number(start+(1-start)*offset)}" ${paint(color,'stop-color','stop-opacity')}/>`);}};
    }
    attributes(stroke){
      const s=this.state;
      let attrs=`transform="matrix(${s.transform.map(number).join(' ')})" opacity="${number(s.globalAlpha)}"`;
      attrs+=stroke?` fill="none" ${paint(s.strokeStyle,'stroke')} stroke-width="${number(s.lineWidth)}" stroke-linejoin="${xml(s.lineJoin)}"${s.dash.length?` stroke-dasharray="${s.dash.map(number).join(' ')}" stroke-dashoffset="${number(s.lineDashOffset)}"`:''}`:` ${paint(s.fillStyle,'fill')} stroke="none"`;
      if(s.shadowBlur>0){
        const id=`shadow-${this.defs.length}`;
        this.defs.push(`<filter id="${id}" x="-100%" y="-100%" width="300%" height="300%" color-interpolation-filters="sRGB"><feGaussianBlur in="SourceAlpha" stdDeviation="${number(s.shadowBlur/2)}" result="blur"/><feOffset in="blur" dx="${number(s.shadowOffsetX)}" dy="${number(s.shadowOffsetY)}" result="offset"/><feFlood ${paint(s.shadowColor,'flood-color','flood-opacity')} result="color"/><feComposite in="color" in2="offset" operator="in" result="shadow"/><feMerge><feMergeNode in="shadow"/><feMergeNode in="SourceGraphic"/></feMerge></filter>`);
        attrs+=` filter="url(#${id})"`;
      }
      return attrs;
    }
    fill(){if(this.path.length)this.elements.push(`<path d="${this.path.join(' ')}" ${this.attributes(false)}/>`);}
    stroke(){if(this.path.length)this.elements.push(`<path d="${this.path.join(' ')}" ${this.attributes(true)}/>`);}
    text(value,x,y,stroke){
      const s=this.state,size=s.font.match(/([\d.]+)px/)?.[1]||'12',weight=s.font.match(/^(\d+)\s/)?.[1]||'400';
      this.elements.push(`<text x="${number(x)}" y="${number(y)}" font-family="system-ui, -apple-system, sans-serif" font-size="${size}" font-weight="${weight}" text-anchor="${s.textAlign==='center'?'middle':s.textAlign==='right'?'end':'start'}" dominant-baseline="${s.textBaseline==='middle'?'central':'alphabetic'}" ${this.attributes(stroke)}>${xml(value)}</text>`);
    }
    fillText(...args){this.text(...args,false);}
    strokeText(...args){this.text(...args,true);}
    mathAxisLabel(order,x,y,size=16,color='#3f4241'){
      if(!window.MajoranaLatex)throw new Error('The LaTeX axis renderer has not loaded.');
      this.elements.push(`<g transform="matrix(${this.state.transform.map(number).join(' ')})" opacity="${number(this.state.globalAlpha)}">${window.MajoranaLatex.toSVG(order,x,y,size,color)}</g>`);
    }
    toSVG(width,height,title,description){
      return `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="${number(width)}" height="${number(height)}" viewBox="0 0 ${number(width)} ${number(height)}" role="img" aria-labelledby="title description"><title id="title">${xml(title)}</title><desc id="description">${xml(description)}</desc><defs>${this.defs.map(d=>typeof d==='function'?d():d).join('')}</defs>${this.elements.join('')}</svg>\n`;
    }
  }
  function download(svg,filename){
    const url=URL.createObjectURL(new Blob([svg],{type:'image/svg+xml;charset=utf-8'}));
    const a=document.createElement('a');a.href=url;a.download=filename;document.body.append(a);a.click();a.remove();
    setTimeout(()=>URL.revokeObjectURL(url),1000);
  }
  window.MajoranaSVG=Object.freeze({Context,download});
})();
