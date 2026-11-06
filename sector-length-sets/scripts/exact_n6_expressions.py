#!/usr/bin/env python3
"""Exact manuscript expressions for the 17 displayed N=6 representatives."""
import cmath
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path

root=Path(__file__).resolve().parents[1]
records={r['id']:r for r in json.loads((root/'data/n6.json').read_text())['states']}
nodes=json.loads((root/'data/n6-manuscript-nodes.json').read_text())
visible={sid for node in nodes for sid in node['stateIds']}
expressions={};vectors={}

def add(sid,latex,amplitudes,*,note='',parameterLatex='',proofReference=''):
    expressions[sid]={'latex':latex,'note':note or 'Dₖ denotes the normalized six-qubit Dicke state |D₆⁽ᵏ⁾⟩.',
                      'proofReference':proofReference}
    if parameterLatex: expressions[sid]['parameterLatex']=parameterLatex
    vector=[0j]*7
    for k,c in amplitudes.items():vector[k]=complex(c)
    vectors[sid]=vector

sqrt=math.sqrt
add('n6-product',r'|0\rangle^{\otimes6}=D_0',{0:1})
add('n6-ghz',r'\displaystyle\frac{D_0+D_6}{\sqrt2}',{0:1/sqrt(2),6:1/sqrt(2)})
add('n6-w',r'D_1',{1:1})
add('n6-octa',r'\displaystyle\frac{D_1+D_5}{\sqrt2}',{1:1/sqrt(2),5:1/sqrt(2)})
add('n6-dicke-2',r'D_2',{2:1})
add('n6-c4',r'\displaystyle\frac{D_0+\sqrt3\,D_4}{2}',{0:.5,4:sqrt(3)/2},proofReference='Spectral equality classification; additional C4 orbit')
add('n6-omega4',r'\displaystyle\frac{\sqrt3\,D_0+\sqrt5\,D_4}{\sqrt8}',{0:sqrt(3/8),4:sqrt(5/8)},proofReference='eq:Omega4-state')
add('n6-psi6',r'\displaystyle\frac{\sqrt{20}\,D_2+D_5}{\sqrt{21}}',{2:sqrt(20/21),5:1/sqrt(21)},
    note='Shares its sector tuple with χ₆, while representing a different collective-rotation orbit.',proofReference='eq:shared-tuple')
add('n6-chi6',r'\displaystyle\frac{2D_0+\sqrt3\,D_5}{\sqrt7}',{0:2/sqrt(7),5:sqrt(3/7)},
    note='Shares its sector tuple with ψ₆, while representing a different collective-rotation orbit.',proofReference='eq:shared-tuple')
add('n6-chi7-10',r'\displaystyle\frac{\sqrt7\,D_0+\sqrt3\,D_5}{\sqrt{10}}',{0:sqrt(.7),5:sqrt(.3)},proofReference='Attained C5 rim, p=7/10')
add('n6-chi1-2',r'\displaystyle\frac{D_0+D_5}{\sqrt2}',{0:1/sqrt(2),5:1/sqrt(2)},proofReference='Attained C5 rim, p=1/2')
add('n6-ac2',r'\displaystyle\frac{\sqrt2\,D_0+\sqrt5\,D_3+\sqrt2\,D_6}{3}',
    {0:sqrt(2)/3,3:sqrt(5)/3,6:sqrt(2)/3},proofReference='eq:C3-family, r=1/9')
add('n6-psi0',r'\displaystyle\frac{D_0+\sqrt2\,D_3+D_6}{2}',
    {0:.5,3:1/sqrt(2),6:.5},proofReference='eq:C3-family, r=0')

u=(sqrt(5)-sqrt(3))/4
add('n6-xi',r'\begin{aligned}\frac12\Bigl[&\sqrt{(1+u_*)(1-2u_*)}\,D_0\\&+\sqrt{(1+u_*)(1+2u_*)}\,D_2\\&+\sqrt{(1-u_*)(1-2u_*)}\,D_4\\&-\sqrt{(1-u_*)(1+2u_*)}\,D_6\Bigr]\end{aligned}',
    {0:sqrt((1+u)*(1-2*u))/2,2:sqrt((1+u)*(1+2*u))/2,4:sqrt((1-u)*(1-2*u))/2,6:-sqrt((1-u)*(1+2*u))/2},
    parameterLatex=r'u_*:=\frac{\sqrt5-\sqrt3}{4}',
    note='Exact normalized representative at the sharp S₁=0, S₃=8 endpoint.',proofReference='eq:Xi')
add('n6-omega',r'\begin{aligned}\frac{1}{4\sqrt{11}}\Bigl[&(\sqrt{21}+3\sqrt2)D_1\\&+7\sqrt2\,D_3\\&+(\sqrt{21}-3\sqrt2)D_5\Bigr]\end{aligned}',
    {1:(sqrt(21)+3*sqrt(2))/(4*sqrt(11)),3:7*sqrt(2)/(4*sqrt(11)),5:(sqrt(21)-3*sqrt(2))/(4*sqrt(11))},
    note='Exact odd-C₂ representative at the endpoint of the sharp-plane equality segment.',proofReference='eq:omega-state')
add('n6-eta6',r'\displaystyle\frac{5D_0+\sqrt{35}\,D_2+\sqrt{15}\,D_4-\sqrt{21}\,D_6}{\sqrt{96}}',
    {0:5/sqrt(96),2:sqrt(35/96),4:sqrt(15/96),6:-sqrt(21/96)},
    note='Sharp paired-null endpoint: S₁=2/27, S₂=65/9, S₃=220/27.',proofReference='Sharp paired-null threshold and equality classification')

# The manuscript's exact cofactor formula evaluated independently at 90 digits.
with localcontext() as ctx:
    ctx.prec=90
    D=Decimal;ell=D(21).sqrt()
    q=D(14027)/D(5375)-D(10647)*ell/D(21500);r=q.sqrt()
    f=((-D(187677)+D(44793)*ell)/D(21500)).sqrt()
    z=(D(49)-D(21)*ell)/D(20);g=D(21)*(D(7)+ell)*r/D(26)
    ap=D(21)/2-D(15)*z/2+D(21)*f/2;am=D(21)/2-D(15)*z/2-D(21)*f/2
    bp=D(21)/2+D(9)*z/2+D(7)*f/2;bm=D(21)/2+D(9)*z/2-D(7)*f/2
    s=D(378)*r/D(13);v=-D(15).sqrt()*g/2;w=s-D(3)*g
    ks=[bp*(bm*am-v*v)-w*w*am,-v*(bm*am-v*v+w*s),v*(w*am+bp*s),s*(w*w-bp*bm)-w*v*v]
    norm=sum(k*k for k in ks).sqrt()
    amplitudes={j:float(-k/norm) for j,k in zip((0,2,4,6),ks)}
add('n6-estar',r'\displaystyle-\frac{k_0D_0+k_2D_2+k_4D_4+k_6D_6}{\sqrt{k_0^2+k_2^2+k_4^2+k_6^2}}',amplitudes,
    parameterLatex=r'''\begin{aligned}
\ell&=\sqrt{21},\quad q_*={14027\over5375}-{10647\ell\over21500},\quad r_*=\sqrt{q_*},\\
f_*&=\sqrt{{-187677+44793\ell\over21500}},\quad \zeta_*={49-21\ell\over20},\\
g_*&={21(7+\ell)\over26}r_*,\\
\alpha_\pm&={21-15\zeta_*\pm21f_*\over2},\quad\beta_\pm={21+9\zeta_*\pm7f_*\over2},\\
s_*&={378\over13}r_*,\quad v_*=-{\sqrt{15}\over2}g_*,\quad w_*=s_*-3g_*,\\
k_0&=\beta_+(\beta_-\alpha_--v_*^2)-w_*^2\alpha_-,\\
k_2&=-v_*(\beta_-\alpha_--v_*^2+w_*s_*),\\
k_4&=v_*(w_*\alpha_-+\beta_+s_*),\\
k_6&=s_*(w_*^2-\beta_+\beta_-)-w_*v_*^2.
\end{aligned}''',
    note='Exact normalized representative of the exceptional aligned C₂ orbit. The auxiliary ζ* is the manuscript’s Qzz parameter, not S₃. The overall minus sign matches the stored state.',
    proofReference='eq:exceptional-C2-state; Proposition prop:exceptional-plane-C2')

assert set(expressions)==visible, (set(expressions)-visible,visible-set(expressions))
audit={}
for sid,vector in vectors.items():
    c=[complex(*v) for v in records[sid]['coefficients']]
    overlap=sum(a.conjugate()*b for a,b in zip(vector,c))
    phase=overlap/abs(overlap)
    error=max(abs(a*phase-b) for a,b in zip(vector,c))
    norm_error=abs(sum(abs(a)**2 for a in vector)-1)
    assert error<1e-9,(sid,error)
    assert norm_error<2e-15,(sid,norm_error)
    audit[sid]={'coefficient_error_up_to_global_phase':error,'normalization_error':norm_error}
text=json.dumps(expressions,ensure_ascii=False,separators=(',',':'))
(root/'data/exact-n6.js').write_text('globalThis.EXACT_N6_EXPRESSIONS = '+text+';\n')
(root/'data/exact-n6.json').write_text(json.dumps(expressions,ensure_ascii=False,indent=2)+'\n')
(root/'data/exact-n6-audit.json').write_text(json.dumps({'status':'PASS','states':len(expressions),'entries':audit},indent=2)+'\n')
print('PASS: 17 exact normalized N6 expressions match the 17 displayed representatives. Max coefficient error:',max(a['coefficient_error_up_to_global_phase'] for a in audit.values()))
