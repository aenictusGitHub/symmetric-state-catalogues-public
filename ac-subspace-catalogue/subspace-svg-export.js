/* Self-contained vector exports of the whole-subspace viewer.
 * Plot geometry remains SVG. KaTeX labels use XHTML with embedded local fonts,
 * following the sector-length atlas export convention; no raster is introduced.
 * compose(viewSvg, {title, subtitle, formulaRows, colorScale, symmetryRows, notes})
 * mathLabel(tex, x, baselineY, fontSize, color, 'start'|'middle'|'end')
 */
(function(global){
  'use strict';
  const NS='http://www.w3.org/2000/svg', XHTML='http://www.w3.org/1999/xhtml';
  const escape=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const number=x=>Number(Number(x).toFixed(3));
  const safeColor=x=>/^(#[\da-f]{3,8}|rgb\([\d.,%\s]+\)|[a-z]+)$/i.test(String(x))?String(x):'#3f4241';
  let sequence=0;
  function math(tex){
    if(!global.katex?.renderToString)throw new Error('LaTeX rendering has not loaded.');
    return global.katex.renderToString(String(tex),{throwOnError:false,strict:'ignore',trust:false,output:'html'});
  }
  function block(html,width,size=14){
    const element=document.createElementNS(XHTML,'div');
    element.setAttribute('style',`box-sizing:border-box;width:${width}px;padding:0;margin:0;color:#343c39;font:${size}px/1.5 "Latin Modern Roman",Georgia,serif;overflow-wrap:break-word;`);
    element.innerHTML=html;
    const host=document.createElement('div');
    host.style.cssText='position:fixed;left:-100000px;top:0;visibility:hidden;pointer-events:none;';
    host.append(element);document.body.append(host);
    // A very long exact expression may need a small reduction to keep every
    // glyph inside the legend rather than clipping at the export boundary.
    for(const label of element.querySelectorAll('[data-svg-math]')){
      const available=width,actual=label.scrollWidth||label.getBoundingClientRect().width;
      if(actual>available)label.style.fontSize=Math.max(.55,available/actual)*size+'px';
    }
    const height=Math.max(1,element.getBoundingClientRect().height);
    host.remove();
    return {element,height};
  }
  const serialized=element=>new XMLSerializer().serializeToString(element);
  const foreign=(element,x,y,width,height)=>`<foreignObject x="${number(x)}" y="${number(y)}" width="${number(width)}" height="${number(height)}">${serialized(element)}</foreignObject>`;
  const mathHtml=tex=>`<span data-svg-math="true" data-tex="${escape(tex)}" style="display:inline-block;max-width:100%;">${math(tex)}</span>`;
  function mathLabel(tex,x,y,size=14,color='#3f4241',anchor='start'){
    const element=document.createElementNS(XHTML,'div');
    element.setAttribute('style',`display:inline-block;white-space:nowrap;margin:0;padding:0;color:${safeColor(color)};font:${size}px/1.2 "Latin Modern Roman",Georgia,serif;`);
    element.innerHTML=mathHtml(tex)+'<span data-svg-baseline="true" style="display:inline-block;width:0;height:0;vertical-align:baseline;"></span>';
    element.querySelector('.katex').style.fontSize='1em';
    const host=document.createElement('div');host.style.cssText='position:fixed;left:-100000px;top:0;visibility:hidden;pointer-events:none;white-space:nowrap;';host.append(element);document.body.append(host);
    const b=element.getBoundingClientRect(),baseline=element.querySelector('[data-svg-baseline]').getBoundingClientRect().top-b.top;
    const width=Math.max(1,b.width)+4,height=Math.max(1,b.height)+4;
    host.remove();
    const left=x-(anchor==='middle'?width/2:anchor==='end'?width:0);
    return foreign(element,left-2,y-baseline-2,width,height);
  }
  function stopsFor(stops){
    return (stops||[]).map(stop=>{
      const offset=Array.isArray(stop)?stop[0]:stop.t??stop.offset;
      const raw=Array.isArray(stop)?stop[1]:stop.color??(stop.rgb?'rgb('+stop.rgb.join(',')+')':'#808080');
      return {offset:Math.max(0,Math.min(1,Number(offset))),color:safeColor(raw)};
    }).filter(stop=>Number.isFinite(stop.offset)).sort((a,b)=>a.offset-b.offset);
  }
  function compose(viewSvg,options={}){
    if(typeof global.ATLAS_EXPORT_FONT_CSS!=='string')throw new Error('Embedded SVG export fonts have not loaded.');
    const doc=new DOMParser().parseFromString(viewSvg,'image/svg+xml'),plot=doc.documentElement;
    if(plot.localName!=='svg'||doc.querySelector('parsererror'))throw new Error('The current view is not valid SVG.');
    if(plot.querySelector('image'))throw new Error('A vector view is required for SVG export.');
    const vb=(plot.getAttribute('viewBox')||'').trim().split(/[ ,]+/).map(Number);
    const width=vb.length===4&&vb.every(Number.isFinite)?vb[2]:parseFloat(plot.getAttribute('width'));
    const height=vb.length===4&&vb.every(Number.isFinite)?vb[3]:parseFloat(plot.getAttribute('height'));
    if(!(width>0&&height>0))throw new Error('The view dimensions are unavailable.');
    const panelWidth=300,panelX=width+40,outputWidth=panelX+panelWidth+26;
    const parts=[];let y=25;
    function add(html,size=14,gap=12,x=panelX,w=panelWidth){
      const b=block(html,w,size);parts.push(foreign(b.element,x,y,w,b.height+3));y+=b.height+gap;
    }
    function heading(text){add(`<div style="font-size:11px;letter-spacing:1.5px;color:#645080;">${escape(text.toUpperCase())}</div>`,11,9);}
    function rows(items){for(const row of items||[]){
      const inset=row.color?25:0,b=block(mathHtml(row.tex)+(row.note?`<div style="font-size:12px;color:#68716b;margin-top:3px;">${escape(row.note)}</div>`:''),panelWidth-inset,14);
      if(row.color)parts.push(`<path d="M${number(panelX)} ${number(y+11)}h16" fill="none" stroke="${safeColor(row.color)}" stroke-width="2.4"/>`);
      parts.push(foreign(b.element,panelX+inset,y,panelWidth-inset,b.height+3));y+=b.height+10;
    }}
    add(`<div style="font-size:19px;line-height:1.3;">${escape(options.title||'Whole-subspace visualization')}</div>`,19,8);
    if(options.subtitle)add(`<div style="color:#68716b;">${escape(options.subtitle)}</div>`,13,15);
    if(options.contextTex)add(mathHtml(options.contextTex),14,16);
    if(options.formulaRows?.length){heading('Multipole weights');rows(options.formulaRows);y+=5;}
    if(options.colorScale){
      const scale=options.colorScale,id='sv-export-scale-'+(++sequence),stops=stopsFor(scale.stops);
      heading(scale.title||'Husimi color scale');
      if(scale.label)add(mathHtml(scale.label),14,7);
      parts.push(`<defs><linearGradient id="${id}" x1="0" y1="0" x2="1" y2="0">${stops.map(stop=>`<stop offset="${stop.offset}" stop-color="${stop.color}"/>`).join('')}</linearGradient></defs><rect x="${panelX}" y="${number(y)}" width="${panelWidth}" height="13" rx="2" fill="url(#${id})" stroke="#b7b7ad" stroke-width=".6"/><text x="${panelX}" y="${number(y+32)}" font-size="12" fill="#343c39" font-family="Latin Modern Roman,Georgia,serif">${escape(scale.min)}</text><text x="${panelX+panelWidth}" y="${number(y+32)}" text-anchor="end" font-size="12" fill="#343c39" font-family="Latin Modern Roman,Georgia,serif">${escape(scale.max)}</text>`);
      if(scale.mid!=null)parts.push(`<text x="${panelX+panelWidth/2}" y="${number(y+32)}" text-anchor="middle" font-size="12" fill="#343c39" font-family="Latin Modern Roman,Georgia,serif">${escape(scale.mid)}</text>`);
      y+=43;if(scale.note)add(`<div style="color:#68716b;">${escape(scale.note)}</div>`,12,14);
    }
    if(options.symmetryRows?.length){heading('Symmetry axes');rows(options.symmetryRows);}
    if(options.notes?.length){y+=3;for(const note of options.notes)add(`<div style="color:#68716b;">${escape(note)}</div>`,12,8);}
    const outputHeight=Math.max(height,y+20);
    plot.setAttribute('x','0');plot.setAttribute('y','0');plot.setAttribute('width',number(width));plot.setAttribute('height',number(height));plot.setAttribute('overflow','hidden');plot.setAttribute('data-layer','view');
    const fontCss=global.ATLAS_EXPORT_FONT_CSS.replace(/\]\]>/g,']]]]><![CDATA[>');
    return `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="${NS}" width="${number(outputWidth)}" height="${number(outputHeight)}" viewBox="0 0 ${number(outputWidth)} ${number(outputHeight)}" role="img"><title>${escape(options.title||'Whole-subspace visualization')}</title><desc>Current vector view with analytical multipole weights where verified, the displayed color scale and symmetry-axis legend. Embedded LaTeX fonts; no raster images.</desc><defs><style><![CDATA[${fontCss}]]></style></defs><rect width="100%" height="100%" fill="#faf8f1"/>${serialized(plot)}<path d="M${width+17} 24V${number(outputHeight-24)}" stroke="#d7d8ce" stroke-width="1"/><g data-layer="legend">${parts.join('')}</g></svg>`;
  }
  global.SubspaceSvgExport={compose,mathLabel};
})(typeof window==='object'?window:globalThis);
