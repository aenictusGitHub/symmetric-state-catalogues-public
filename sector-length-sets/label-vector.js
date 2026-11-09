/* Convert the visual KaTeX layout to native SVG outlines. */
(function (global) {
  'use strict';
  const escape = value => String(value).replace(/[&<>"']/g, c =>
    ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&apos;'})[c]);
  const number = value => Number(value.toFixed(6));
  function glyph(fonts, style, character) {
    const italic = style.fontStyle !== 'normal', bold = Number(style.fontWeight) >= 600;
    const variant = bold ? (italic ? 'BoldItalic' : 'Bold') : (italic ? 'Italic' : 'Regular');
    const families = style.fontFamily.split(',').map(f => f.trim().replace(/["']/g,''));
    for (const family of [...families, 'KaTeX_Main', 'Latin Modern Roman']) {
      for (const suffix of [variant,'Regular']) {
        const font = fonts[family+'-'+suffix];
        if (font && Object.hasOwn(font[1],character)) return [font[0],font[1][character]];
      }
    }
    return null;
  }

  function toSVG(element,x,y) {
    const fonts=global.ATLAS_LABEL_OUTLINES;
    if(!fonts)return null;
    const clone=element.cloneNode(true), style=getComputedStyle(element);
    clone.removeAttribute('id');
    clone.querySelectorAll('.katex-mathml').forEach(node=>node.remove());
    Object.assign(clone.style,{position:'fixed',left:'-10000px',top:'0',transform:'none',
      width:'max-content',maxWidth:'none',pointerEvents:'none',color:style.color,
      fontFamily:style.fontFamily,fontSize:style.fontSize,lineHeight:style.lineHeight,
      fontStyle:style.fontStyle,fontWeight:style.fontWeight});
    element.parentNode.append(clone);
    try {
      const box=clone.getBoundingClientRect(),paths=[],walker=document.createTreeWalker(clone,NodeFilter.SHOW_TEXT);
      const range=document.createRange(),nodes=[];
      for(let node=walker.nextNode();node;node=walker.nextNode())nodes.push(node);
      for(const node of nodes){
        if(!node.data.trim())continue;
        const parent=node.parentElement, font=getComputedStyle(parent);
        const marker=document.createElement('span');
        marker.style.cssText='display:inline-block;width:0;height:0;padding:0;margin:0;vertical-align:baseline';
        node.after(marker);
        const baseline=marker.getBoundingClientRect().top-box.top;
        marker.remove();
        let index=0;
        for(const c of node.data){
          const start=index;index+=c.length;
          if(/[\s\u200B]/u.test(c))continue;
          const outline=glyph(fonts,font,c);
          if(!outline)return null;
          if(!outline[1])continue;
          range.setStart(node,start);range.setEnd(node,index);
          const rect=range.getBoundingClientRect(),scale=parseFloat(font.fontSize)/outline[0];
          paths.push('<path d="'+escape(outline[1])+'" fill="'+escape(font.color)+'" transform="translate('+
            number(rect.left-box.left)+' '+number(baseline)+') scale('+number(scale)+' '+number(-scale)+')"/>');
        }
      }
      // KaTeX uses SVG for stretchy delimiters and CSS borders for fraction rules.
      for(const graphic of clone.querySelectorAll('.katex-html svg')){
        const rect=graphic.getBoundingClientRect(),svg=graphic.cloneNode(true);
        for(const [key,val] of Object.entries({x:rect.left-box.left,y:rect.top-box.top,width:rect.width,height:rect.height})){
          svg.setAttribute(key,number(val));svg.style.removeProperty(key);
        }
        svg.setAttribute('color',getComputedStyle(graphic).color);
        paths.push(new XMLSerializer().serializeToString(svg));
      }
      for(const rule of clone.querySelectorAll('.frac-line,.sqrt-line,.overline-line,.underline-line')){
        const rect=rule.getBoundingClientRect(),ruleStyle=getComputedStyle(rule);
        const height=parseFloat(ruleStyle.borderBottomWidth)||parseFloat(ruleStyle.borderTopWidth);
        if(height>0)paths.push('<rect x="'+number(rect.left-box.left)+'" y="'+number(rect.top-box.top)+
          '" width="'+number(rect.width)+'" height="'+number(height)+'" fill="'+escape(ruleStyle.color)+'"/>');
      }
      if(!paths.length)return null;
      const source=[...element.querySelectorAll('annotation[encoding="application/x-tex"]')].map(node=>node.textContent).join(', ');
      return '<g data-tex="'+escape(source||clone.textContent)+'" aria-label="'+escape(source||clone.textContent)+
        '" opacity="'+escape(style.opacity)+'" transform="translate('+number(x-box.width/2)+' '+
        number(y-box.height/2)+')">'+paths.join('')+'</g>';
    }finally{clone.remove();}
  }
  global.AtlasLabelVector=Object.freeze({toSVG});
})(globalThis);
