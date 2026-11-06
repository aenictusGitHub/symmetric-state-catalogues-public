/* Inline mathematics is authored in $...$ and typeset locally with KaTeX. */
(function(){
  const stateTex={
    'n3-product':String.raw`|0\rangle^{\otimes3}`,'n3-ghz':String.raw`|\mathrm{GHZ}(3)\rangle`,'n3-w':String.raw`|W(3)\rangle`,
    'n4-product':String.raw`|0\rangle^{\otimes4}`,'n4-ghz':String.raw`|\mathrm{GHZ}(4)\rangle`,'n4-tetra':String.raw`|\mathrm{tetra}\rangle`,'n4-w':String.raw`|W(4)\rangle`,'n4-dicke2':String.raw`|D_4^{(2)}\rangle`,
    'n5-product':String.raw`|0\rangle^{\otimes5}`,'n5-ghz':String.raw`|\mathrm{GHZ}(5)\rangle`,'n5-w':String.raw`|W(5)\rangle`,'n5-dicke2':String.raw`|D_5^{(2)}\rangle`,'n5-psia':String.raw`|\psi_5^{(a)}\rangle`,'n5-psib':String.raw`|\psi_5^{(b)}\rangle`,'n5-psic':String.raw`|\psi_5^{(c)}\rangle`,
    'n6-product':String.raw`|0\rangle^{\otimes6}`,'n6-ghz':String.raw`|\mathrm{GHZ}(6)\rangle`,'n6-w':String.raw`|W(6)\rangle`,'n6-octa':String.raw`|\mathrm{octa}\rangle`,'n6-dicke-2':String.raw`|D_6^{(2)}\rangle`,'n6-dicke-3':String.raw`|D_6^{(3)}\rangle`,'n6-c4':String.raw`|C_4\rangle`,'n6-omega4':String.raw`|\Omega_4\rangle`,'n6-psi6':String.raw`|\psi_6\rangle`,'n6-chi6':String.raw`|\chi_6\rangle`,'n6-chi7-10':String.raw`|\chi_{7/10}\rangle`,'n6-chi1-2':String.raw`|\chi_{1/2}\rangle`,'n6-ac2':String.raw`|\mathrm{AC}_2\rangle`,'n6-psi0':String.raw`|\psi_0\rangle`,'n6-xi':String.raw`|\Xi_*\rangle`,'n6-omega':String.raw`|\omega\rangle`,'n6-estar':String.raw`|E_*\rangle`,'n6-eta6':String.raw`|\eta_6\rangle`
  };
  const escape=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  function html(text){return String(text).split(/(\$[^$]+\$)/g).map(s=>s.startsWith('$')&&s.endsWith('$')?katex.renderToString(s.slice(1,-1),{throwOnError:false,strict:false}):escape(s)).join('');}
  function state(id,fallback){if(/^n5-top-/.test(id||''))return '$'+String.raw`|\mathrm{GHZ}(5,\phi)\rangle`+'$';return stateTex[id]?'$'+stateTex[id]+'$':String(fallback||'').replace(/C([₂₃₄₅₆])/g,(_,s)=>'$C_'+({'₂':2,'₃':3,'₄':4,'₅':5,'₆':6}[s])+'$');}
  function axis(text){const digits={'₁':1,'₂':2,'₃':3,'₄':4,'₅':5,'₆':6};return /^S[₁₂₃₄₅₆]$/.test(text)?'$S_'+digits[text[1]]+'$':text;}
  function renderStatic(){if(window.renderMathInElement)renderMathInElement(document.body,{delimiters:[{left:'$$',right:'$$',display:true},{left:'$',right:'$',display:false}],throwOnError:false,strict:false});}
  window.AtlasMath={html,state,axis,renderStatic};
})();
