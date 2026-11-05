#!/usr/bin/env python3
"""Audit the finite dimension table and present proved patterns separately from conjectures."""
from pathlib import Path
import html, json, re, hashlib
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def m(x,display=True):
    tag='div' if display else 'span'
    return f'<{tag} class="math" data-display="{str(display).lower()}" data-tex="{html.escape(x,quote=True)}"></{tag}>'
def p(x):return '<p>'+x+'</p>'
def link(url,label):return f'<a href="{html.escape(url,quote=True)}">{label} ↗</a>'
def refs(items):return '<div class="proof-links">'+''.join(link(a,b) for a,b in items)+'</div>'

def main():
    raw=(ROOT/'data/catalogue.json').read_bytes();data=json.loads(raw)
    rows={(e['n'],e['t']):e for e in data['entries']}
    exact=[e for e in data['entries'] if e['lower']==e['upper']]
    positive=[e for e in exact if e['lower']>0]
    empty=[e for e in exact if e['lower']==0]
    assert rows[16,6]['upper']==0
    plateau=[]
    for n in data['numbers']:
        run=[]
        for t in range(2,9):
            e=rows.get((n,t))
            if e and e['lower']==e['upper'] and e['lower']>0:
                if run and e['lower']!=rows[n,run[-1]]['lower']:
                    if len(run)>1:plateau.append({'n':n,'orders':run,'k':rows[n,run[-1]]['lower']})
                    run=[]
                run.append(t)
            else:
                if len(run)>1:plateau.append({'n':n,'orders':run,'k':rows[n,run[-1]]['lower']})
                run=[]
    third=[e for e in positive if e['t']==3 and e['n']>=10]
    third_matches=[e for e in third if e['lower']==(e['n']-1)//6]
    third_exceptions=[{'n':e['n'],'provedMaximum':e['lower'],'fittedValue':(e['n']-1)//6} for e in third if e['lower']!=(e['n']-1)//6]
    assert (55+7)//7 < (55-1)//6
    findings={'catalogueSHA256':hashlib.sha256(raw).hexdigest(),'entries':len(rows),
        'positiveExactMaxima':len(positive),'provedEmpty':len(empty),'unsettled':len(rows)-len(exact),
        'exactPositiveCountByOrder':{str(t):sum(e['t']==t for e in positive) for t in data['orders']},
        'plateausAtOrdersAtLeastTwo':plateau,
        'thirdOrderFiniteFit':{'formula':'floor((N-1)/6)','minimumN':10,'matchingProvedSizes':[e['n'] for e in third_matches],
            'exactCounterexamples':third_exceptions,
            'universalExtensionRuledOutAtN':55,'fittedValue':9,'newRigorousUpper':8},
        'newGeneralBounds':[{'minimumN':6,'minimumOrder':3,'bound':'7k<=N+7'},
            {'minimumN':9,'minimumOrder':4,'bound':'10k<=N+10'},
            {'minimumN':16,'minimumOrder':4,'bound':'7k<=N-5'},
            {'minimumN':19,'minimumOrder':3,'bound':'6k<=N-4'}],
        'additionalSettledCells':[{'n':18,'t':4,'previousInterval':[1,2],'newInterval':[1,1]}],
        'newSettledCell':{'n':16,'t':6,'previousInterval':[0,1],'newInterval':[0,0]},
        'entriesData':[{'n':e['n'],'t':e['t'],'lower':e['lower'],'upper':e['upper'],'status':e['status']} for e in data['entries']]}
    tight=json.loads((HERE/'packing-tightness-verification.json').read_text());assert tight['passed']
    findings['unroundedPackingEquality']=tight['classification']
    assert tight['classification']=={'t1':[[2,1],[4,2],[6,3]],'t2':[[7,2]],'t3':[[6,1]],'t4':[]}
    assert all(e['upper']<=(e['n']-1)//6 for e in data['entries'] if e['t']>=3 and e['n']>=7)
    assert all(e['upper']<=e['n']//9 for e in data['entries'] if e['t']>=4)
    (HERE/'patterns.json').write_text(json.dumps(findings,indent=2)+'\n')
    css=re.search(r'<style>(.*?)</style>',(ROOT/'research/t2-qec-upper/index.html').read_text(),re.S).group(1)
    css+='''
    .math{max-width:100%;overflow-x:auto}.matrix{min-width:620px}.matrix td,.matrix th{text-align:center;padding:.65rem .8rem}.matrix td.exact{background:#e4f0e7;color:#185e50;font-weight:650}.matrix td.zero{background:#efebf5;color:#695380}.matrix td.open{color:#776d60}.matrix a{color:inherit;text-decoration:none}.badge{display:inline-block;font-size:.78rem;letter-spacing:.04em;text-transform:uppercase;padding:.25rem .55rem;border:1px solid #a7c3b4;border-radius:4px;margin-bottom:.6rem;color:#185e50}.caution{border-color:#cbb788;color:#7b5e1a}.scroll{overflow-x:auto}.legend{font-size:.9rem;color:var(--muted,#626b67)}.equation-pair{display:grid;grid-template-columns:1fr 1fr;gap:1rem}.equation-pair>div{background:#e9eee9;padding:1rem;border:1px solid #cdd5cd;border-radius:4px}summary{cursor:pointer;font-weight:600}details{margin:1.3rem 0}li{margin:.5rem 0}section{scroll-margin-top:1rem}@media(max-width:660px){.equation-pair{grid-template-columns:1fr}.katex{font-size:1.03em}}'''
    b='<p class="eyebrow">Anticoherent subspaces · Exact maxima across orders and sizes</p>'
    b+=link('../../index.html','← AC-subspace catalogue')+'<h1>Patterns that survive exact verification</h1>'
    b+=p(f'The complete catalogue contains {len(rows)} dimension intervals: {len(positive)} positive exact maxima, {len(empty)} proved empty cases and {len(rows)-len(exact)} unresolved intervals. Comparing all of them reveals algebraic links and equality obstructions, as well as a tempting numerical fit that cannot continue indefinitely.')
    b+='<div class="result"><span class="badge">New proved results</span><div class="equation-pair"><div>'+m(r'7K(N,3)\le N+7\quad(N\ge6)')+'</div><div>'+m(r'10K(N,4)\le N+10\quad(N\ge9)')+'</div></div>'+m(r'\boxed{K(16,6)=0}')+p('The two general bounds improve the large-size theory. The last result closes the displayed sixteen-qubit sixth-order gap. All three apply to arbitrary complex states or subspaces.')+'</div>'
    b+='<nav class="jump-links" aria-label="On this page"><a href="#bounds">Denominator pattern</a><a href="#bridges">Exact projector links</a><a href="#equality">Equality obstructions</a><a href="#plateaus">Plateaus and exceptions</a><a href="#table">Complete table</a></nav>'

    b+='<section id="bounds"><h2>A new denominator pattern in rigorous upper bounds</h2>'
    b+=p('Let '+m(r'D_t=\lfloor(t+2)^2/4\rfloor',False)+'. The optimized rectangular-tensor counts are '+m(r'D_t=2,4,6,9,12,\ldots',False)+', the shifted '+link('https://oeis.org/A002620','quarter-square sequence A002620')+'. The auxiliary-spin method now improves these first four denominators by one:')
    b+='<div class="scroll"><table><thead><tr><th>AC order</th><th>Tensor count '+m(r'D_t',False)+'</th><th>Auxiliary dimension</th><th>Proved bound</th><th>Domain</th></tr></thead><tbody>'
    for t,d,a,bound,domain in [(1,2,3,r'3k\le N+3',r'N\ge2'),(2,4,5,r'5k\le N+5',r'N\ge4'),(3,6,7,r'7k\le N+7',r'N\ge6'),(4,9,10,r'10k\le N+10',r'N\ge9')]:
        b+=f'<tr><td>{t}</td><td>{d}</td><td>{a}</td><td>'+m(bound,False)+'</td><td>'+m(domain,False)+'</td></tr>'
    b+='</tbody></table></div>'
    b+=p('The last two rows are new. For each order, choose a polynomial of degree at most '+m(r't',False)+' in a physical–auxiliary spin coupling. Whole-subspace AC makes its centered compression zero. A polynomial with only the highest-spin sector positive bounds the code dimension by the dimension of that sector.')
    b+='<details><summary>Third-order polynomial and proof</summary>'+m(r'j=N/2,\quad X=\mathbf J\cdot\mathbf S,\quad S=3,')+m(r'W=(X+2j+5)(X+j+6)(X-2j+3)-15(2j^2+3j+6)I.')
    b+=p('AC3 gives '+m(r'P X P=0,\ P X^2P=4j(j+1)P,\ P X^3P=-2j(j+1)P',False)+' with '+m(r'P=P_{\mathcal C}\otimes I_7',False)+', so '+m(r'PWP=0',False)+'. All sectors below total spin '+m(r'j+3',False)+' are strictly negative, while the highest sector is positive. Projection into the positive sector is injective on the '+m(r'7k',False)+'-dimensional code–auxiliary space. Hence '+m(r'7k\le N+7',False)+'.')
    b+=refs([('spin3-proof.tex','Full proof and spectrum'),('spin3-verification.json','Exact symbolic and physical checks')])+'</details>'
    b+='<details><summary>Fourth-order polynomial and proof</summary>'+m(r'S=9/2,\qquad F(x)=\frac{(2x-7j+9)(2x+j+25)(2x+3j+24)(2x+9j+9)}{16},')+m(r'c_j=\frac{9(j+1)(4j^3+724j^2+2637j+6750)}{20},\quad W=F(X)-c_j I.')
    b+=p('The scalar '+m(r'c_j',False)+' is the normalized physical trace. At every lower spectral node '+m(r'F\le0',False)+', whereas '+m(r'c_j>0',False)+'. The highest eigenvalue of '+m(r'W',False)+' is strictly positive. AC4 implies zero compression on the '+m(r'10k',False)+'-dimensional code–auxiliary space, which must inject into the '+m(r'N+10',False)+'-dimensional highest-spin sector.')
    b+=refs([('spin4-proof.tex','Full proof and polynomial'),('spin4-independent-audit.tex','Independent Cartesian-moment derivation'),('spin4-independent-verification.json','Independent exact audit')])+'</details>'
    b+=p('Combining these with the existing bounds gives, for example, '+m(r'K(55,3)\le8',False)+' instead of nine, and '+m(r'K(99,4)\le10',False)+' instead of eleven. Neither inequality alone changes a displayed cell, since the catalogue currently ends at selected sizes up to 42.')
    b+='<details open><summary>A complementary fourth-order bound closes the eighteen-qubit gap</summary>'+m(r'7K(N,t)\le N-5\qquad(N\ge16,\ t\ge4).')+p('An auxiliary spin-three quartic is positive only on the lowest total-spin sector, of dimension '+m(r'N-5',False)+'. Zero compression of the code–auxiliary space gives this bound. At eighteen qubits, '+m(r'7k\le13',False)+' proves '+m(r'K(18,4)=1',False)+'. At nineteen qubits it still permits dimension two.')+refs([('../n18-t4/index.html','Short proof and exact witness'),('../n18-t4/proof.tex','General analytical theorem')])+'</details>'
    b+='<details open><summary>A cubic bound closes the nineteen- and twenty-qubit third-order gaps</summary>'+m(r'6K(N,t)\le N-4\qquad(N\ge19,\ t\ge3).')+p('An auxiliary spin five-halves gives a cubic with positive spectrum only in the lowest sector, of dimension '+m(r'N-4',False)+'. Consequently third-order dimension cannot exceed two at nineteen or twenty qubits. The explicit analytical planes attain this bound.')+refs([('../t3-spin5half-inertia/index.html','General cubic proof'),('../n19-t3/index.html','Nineteen-qubit plane'),('../n20-t3/index.html','Twenty-qubit plane')])+'</details>'
    b+='<span class="badge caution">Unproved extension</span>'+p('The coefficients '+m(r'3,5,7,10',False)+' equal '+m(r'D_t+1',False)+' for '+m(r't=1,2,3,4',False)+'. This suggests testing '+m(r'(D_t+1)k\le N+D_t+1',False)+' at higher orders. No such general-order theorem or matching construction is asserted.')
    b+='<details><summary>The third-order sequence fit has exact counterexamples</summary>'+p('Among the proved positive maxima at '+m(r'N\ge10',False)+', the expression '+m(r'\lfloor(N-1)/6\rfloor',False)+' agrees at sizes '+', '.join(str(e['n']) for e in third_matches)+'. It fails at '+', '.join('N='+str(e['n'])+' (maximum '+str(e['provedMaximum'])+', fitted value '+str(e['fittedValue'])+')' for e in third_exceptions)+'. These are exact finite counterexamples, established by subspace-exclusion certificates. Independently, the expression predicts nine at 55 qubits, contradicting the upper bound eight. Unresolved cells must not be discarded when assessing a pattern.')+refs([('../n14-t3/index.html','Fourteen-qubit exclusion'),('../n15-t3/index.html','Fifteen-qubit exclusion')])+'</details></section>'

    b+='<section id="bridges"><h2>Three optimal objects linked by an exact projector identity</h2>'
    b+=p('The optimal seven-qubit, two-dimensional AC2 subspace has basis')
    b+=m(r'u=\sqrt{3/10}|D_7^{(0)}\rangle+\sqrt{7/10}|D_7^{(5)}\rangle,\quad v=\sqrt{7/10}|D_7^{(2)}\rangle-\sqrt{3/10}|D_7^{(7)}\rangle.')
    b+=p('Let '+m(r'P_7=|u\rangle\langle u|+|v\rangle\langle v|',False)+'. For the catalogue’s optimal six-qubit, three-dimensional AC1 projector '+m(r'P_6',False)+', exact substitution gives')
    b+=m(r'\boxed{P_6=I_7-2\operatorname{Tr}_1(P_7).}')
    b+=p('Thus the kernel of the reduced normalized code state is an optimal AC1 subspace: '+m(r'\operatorname{Tr}_1(P_7/2)=(I_7-P_6)/4',False)+'. This is a concrete reduced-state construction; the relevant subspace is its kernel.')
    b+=p('The centered projector '+m(r'Y=P_6-3I_7/7',False)+' has only rank-six spin tensors. Time-reversed vectorization gives')
    b+=m(r'\boxed{|\Phi_{12}\rangle=\frac{-\sqrt7|D_{12}^{(1)}\rangle+\sqrt{11}|D_{12}^{(6)}\rangle+\sqrt7|D_{12}^{(11)}\rangle}{5},\qquad\rho_5=I_6/6.}')
    b+=p('The result is the icosahedral state. Hence the exact maxima '+m(r'(N,t,k)=(6,1,3),(7,2,2),(12,5,1)',False)+' share an explicit algebraic link. A second link takes the optimal '+m(r'(4,1,2)',False)+' projector to the cubic eight-qubit AC3 state:')
    b+=m(r'|\Phi_8\rangle=\frac{\sqrt{30}}{12}\bigl(|D_8^{(0)}\rangle+|D_8^{(8)}\rangle\bigr)+\frac{\sqrt{21}}6|D_8^{(4)}\rangle.')
    b+='<details><summary>The general sufficient construction and its limitation</summary>'+p('If a nontrivial projector '+m(r'P',False)+' on '+m(r'n',False)+' symmetric qubits has its centered part '+m(r'Y=P-\operatorname{rank}(P)I/(n+1)',False)+' entirely in tensor rank '+m(r'n',False)+', then its time-reversed vectorization is a symmetric '+m(r'2n',False)+'-qubit pure state. Since '+m(r'Y^2',False)+' is a linear combination of '+m(r'I',False)+' and '+m(r'Y',False)+', its '+m(r'n-1',False)+'-qubit reduction is maximally mixed: the state is '+m(r'(n-1)',False)+'-AC.')
    b+=p('This special projector condition cannot hold for '+m(r'n\ge7',False)+': the resulting pure state would violate the rectangular-tensor bound. These exact links are exceptional finite constructions, not an infinite family and not a classification of optimal codes.')+'</details>'
    b+=refs([('bridge-proof.tex','Construction proof and all bases'),('bridge-verification.json','Exact projectors and cross reductions'),('bridge-verify.py','Reproduce the calculation')])+'</section>'

    b+='<section id="equality"><h2>Tight packing imposes extra multipole equations</h2>'
    b+=p('For pure states, the tensor-counting thresholds have the same numerical form as the '+link('https://www.math.unipd.it/~alvise/PHD_2024/LEZIONI/SPHERE/PDF/DELSARTE%20Spherical%20codes%20and%20designs.pdf','Fisher bounds for spherical designs')+':')
    b+=m(r'N\ge(r+1)^2\quad(t=2r),\qquad N\ge(r+1)(r+2)\quad(t=2r+1).')
    b+=p('At orders one, two, three and five, familiar antipodal, tetrahedral, octahedral and icosahedral representatives attain these thresholds. This is useful structural motivation. Anticoherence and a spherical-design Majorana constellation are not assumed equivalent; '+link('https://arxiv.org/abs/1509.08300','the point-group study')+' distinguishes the spin-state conditions.')
    b+='<div class="result"><span class="badge">Newly settled cell</span>'+m(r'K(16,6)=0')+p('At '+m(r'N=16,t=6',False)+', fifteen rank-one, rank-two and rank-three tensor images would fill the spin-seven output space. Equality therefore forces the channel output to be '+m(r'I_{15}/17',False)+'. Its rank-seven multiplier is')+m(r'\frac{3\sqrt{69}}{170}+\frac{\sqrt{69}}{306}-\frac{7679\sqrt{69}}{377910}=\frac{5\sqrt{69}}{8398}\ne0.')+p('Consequently every rank-seven multipole must vanish as well. Seventh-order AC would require twenty orthogonal images in dimension sixteen, which is impossible. This upgrades the former interval '+m(r'0\le K(16,6)\le1',False)+' to an exact zero.')+'</div>'
    b+=p('The same exact criterion rules out the tight even-order pairs '+m(r'(N,t)=(9,4),(25,8),(36,10),(49,12),(64,14),(81,16)',False)+'. These are finitely verified exclusions. The required multiplier has not been proved nonzero at every order.')
    b+=refs([('even-tight-obstruction.tex','Proof and general conditional criterion'),('even-tight-verification.json','Exact channel coefficients'),('even-tight-independent-verification.json','Independent Clebsch–Gordan audit')])
    b+='<h3>Complete packing-equality classification through order four</h3>'
    b+=p('The new auxiliary bounds restrict equality in the original unrounded packing inequalities to a finite list of sizes. Exact channel multipliers and the existing seventeen-qubit obstruction settle every remaining case:')
    b+='<div class="scroll"><table><thead><tr><th>Order</th><th>Unrounded equality</th><th>All attaining '+m(r'(N,k)',False)+'</th></tr></thead><tbody>'
    for t,equality,cases in [(1,r'2k=N',r'(2,1),(4,2),(6,3)'),(2,r'4k=N+1',r'(7,2)'),(3,r'6k=N',r'(6,1)'),(4,r'9k=N+1',r'\text{none}')]:
        b+=f'<tr><td>{t}</td><td>'+m(equality,False)+'</td><td>'+m(cases,False)+'</td></tr>'
    b+='</tbody></table></div>'+m(r'6K(N,3)\le N-1\quad(N\ge7),\qquad9K(N,4)\le N\quad(N\ge4).')
    b+=p('This classifies equality in the unrounded tensor-counting inequalities. It does not classify saturation after integer rounding or assert attainment of the stronger upper bounds. The strict refinements agree with the already applied bounds in every displayed cell.')
    b+=refs([('packing-tightness.tex','Complete equality classification and proof'),('packing-tightness-verification.json','Exact finite exclusions and attaining witnesses'),('packing-tightness-verify.py','Independent verifier')])
    b+='<details><summary>Third-order saturation leads to a Pell-type condition</summary>'+p('The existing third-order equality criterion says that '+m(r'6k=N\ge12',False)+' requires an integer rank '+m(r'4\le L\le N-1',False)+' satisfying')+m(r'2L(L+1)=(N-1)(N+2)\quad\Longleftrightarrow\quad(2N+1)^2-2(2L+1)^2=7.')+p('This makes equality arithmetically exceptional. A solution is necessary, not sufficient. The new '+m(r'7k\le N+7',False)+' bound actually excludes '+m(r'6k=N',False)+' for every '+m(r'N>42',False)+', even when the Pell-type condition has a solution.')+refs([('../general-bounds/index.html#equality','Earlier equality proof')])+'</details></section>'

    b+='<section id="plateaus"><h2>Plateaus, size exceptions and a limit on nesting</h2>'
    b+=p('Several plateaus have a direct explanation: the same higher-order witness already reaches the lower-order upper bound.')
    b+=m(r'K(12,3)=K(12,4)=K(12,5)=1,\qquad K(13,2)=K(13,3)=2,')+m(r'K(16,4)=K(16,5)=1.')
    b+=p('Increasing the number of qubits need not increase the maximum:')+m(r'K(7,2)=2>K(8,2)=K(9,2)=1.')
    b+=p('In contrast, increasing the order always decreases or preserves the maximum. A smallest known size is therefore not an automatic existence threshold for all larger sizes.')
    b+=p('There is also a rigorous limit on construction by nesting. If an '+m(r's',False)+'-AC parent has dimension '+m(r'd',False)+' and contains a '+m(r'2s',False)+'-AC subspace of dimension '+m(r'k',False)+', then')
    b+=m(r'\boxed{d+s(s+2)k\le N+1.}')
    b+=p('The rank-one through rank-'+m(r's',False)+' tensor images of the smaller code are mutually orthogonal and orthogonal to the entire parent. At seven qubits, '+m(r's=1,k=2',False)+' forces '+m(r'd\le2',False)+', although '+m(r'K(7,1)=3',False)+'. Thus no optimal seven-qubit AC2 code sits inside a maximum-dimensional AC1 parent.')
    b+=refs([('nesting-proof.tex','Nesting inequality and proof')])+'</section>'

    b+='<section id="table"><h2>Every catalogue entry included in the comparison</h2>'+p('Green cells are positive exact maxima; purple zeros are proved nonexistence. An interval remains unresolved. Numerical or conjectured dimensions are not substituted for exact maxima. Select a cell to inspect its proofs and witnesses.')
    b+='<div class="scroll"><table class="matrix"><thead><tr><th scope="col">'+m(r'N\backslash t',False)+'</th>'+''.join(f'<th scope="col">{t}</th>' for t in data['orders'])+'</tr></thead><tbody>'
    for n in data['numbers']:
        b+=f'<tr><th scope="row">{n}</th>'
        for t in data['orders']:
            e=rows.get((n,t))
            if not e:b+='<td>—</td>';continue
            cls='zero' if e['upper']==0 else 'exact' if e['lower']==e['upper'] else 'open'
            text=str(e['lower']) if e['lower']==e['upper'] else f"{e['lower']}–{e['upper']}"
            b+=f'<td class="{cls}"><a href="../../index.html#n={n}&amp;t={t}" aria-label="N {n}, order {t}, dimension {text}">{text}</a></td>'
        b+='</tr>'
    b+='</tbody></table></div>'+refs([('patterns.json','Complete machine-readable pattern audit'),('../../data/catalogue.json','Full catalogue data')])+'</section>'
    b+='<footer>'+link('../../index.html','Return to the AC-subspace catalogue')+'<span>Exact identities · rigorous bounds · open extensions labelled separately</span></footer>'
    head='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Patterns across exact AC-subspace maxima</title><link rel="stylesheet" href="../../vendor/katex/katex.min.css"><script defer src="../../vendor/katex/katex.min.js"></script><style>'+css+'</style><script>document.addEventListener("DOMContentLoaded",()=>document.querySelectorAll("[data-tex]").forEach(el=>{if(window.katex)el.innerHTML=katex.renderToString(el.dataset.tex,{displayMode:el.dataset.display==="true",throwOnError:false,strict:"ignore",trust:false});}));</script></head><body><main>'
    (HERE/'index.html').write_text(head+b+'</main></body></html>\n')
    print(json.dumps({k:v for k,v in findings.items() if k!='entriesData'},indent=2))

if __name__=='__main__':main()
