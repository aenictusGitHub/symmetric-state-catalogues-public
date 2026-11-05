"""Build the offline Majorana-degeneracy catalogue from exact research records."""
from pathlib import Path
from collections import Counter
import json, math, cmath, csv, re, html
import sympy as s
from partition_catalogue import reindex
from liege_library import attach_library
from acliege_exact import attach_acliege_exact

ROOT=Path(__file__).resolve().parents[1]
def read(path):return json.loads((ROOT/path).read_text())
def body(text):
    return re.sub(r'\\\((.*?)\\\)',lambda m:'<span class="math" data-tex="'+html.escape(m.group(1),quote=True)+'" data-display="false">'+html.escape(m.group(1))+'</span>',text)
def formula(n,coeff):
    parts=[]
    for r,value in enumerate(coeff):
        c=s.sympify(value)
        if c==0:continue
        negative=c.could_extract_minus_sign()
        if negative:c=-c
        a='' if c==1 else s.latex(c)
        if c.is_Add:a=r'\left('+a+r'\right)'
        parts.append((' - ' if negative else ' + ')+a+rf'|D_{{{n}}}^{{({r})}}\rangle')
    return ''.join(parts).removeprefix(' + ').removeprefix(' - ') if not parts[0].startswith(' - ') else '-'+''.join(parts)[3:]
def point(z,m):
    # P(z)=prod(a_i-b_i*z), so a physical spinor has b_i/a_i=1/z.
    if z is None:return {'x':0.,'y':0.,'z':1.,'multiplicity':m}
    q=1+abs(z)**2
    return {'x':2*z.real/q,'y':-2*z.imag/q,'z':(abs(z)**2-1)/q,'multiplicity':m}
def roots_for(record):
    result=[]
    if record['zeroMultiplicity']:result.append(point(0j,record['zeroMultiplicity']))
    if record['infinityMultiplicity']:result.append(point(None,record['infinityMultiplicity']))
    step=record['nonzeroRootSubstitution']['step'];y=s.Symbol('y')
    for factor in record['squarefreeFactors']:
        p=s.Poly(s.sympify(factor['polynomial']),y)
        roots=s.nroots(p,n=30,maxsteps=500)
        for root in roots:
            value=complex(root);angle=cmath.phase(value);radius=abs(value)**(1/step)
            for k in range(step):result.append(point(cmath.rect(radius,(angle+2*math.pi*k)/step),factor['multiplicity']))
    assert sum(p['multiplicity'] for p in result)==record['n']
    return result

ac=read('research/ac-catalogue.json')
multiplicity={r['stateId']:r for r in read('research/multiplicities.json')['states']}
states=[]
for raw in read('research/states-input.json')['states']:
    st=dict(raw);r=multiplicity[st['id']]
    st.update(partition=r['partition'],maxMultiplicity=r['maxMultiplicity'],distinctPoints=r['distinctPoints'],
              polynomialLatex=st.get('factoredPolynomialLatex',r['homogeneousPolynomialLatex']),majoranaPoints=roots_for(r),stateLatex=formula(st['n'],st['exactCoefficients']))
    st['coefficientLatex']=[s.latex(s.sympify(c)) for c in st['exactCoefficients']]
    st['verification']={'method':'Exact symbolic reduced matrices and algebraic squarefree factorization, including the projective point at infinity.','exactAC':True,'exactMultiplicity':True,'certificatePath':'./research/multiplicities.json','acCertificatePath':'./research/states-verification.json'}
    st['notes']=st.get('notes','')
    states.append(st)

# This state is defined exactly by a nonsingular rational system in a certified box.
w=read('research/n11-g2/witness.json');center=list(map(s.Rational,w['allCoordinateCenter']))
coeff=[[float(s.sqrt(s.binomial(11,r))*center[r]),float(s.sqrt(s.binomial(11,r))*center[r+10])] for r in range(10)]+[[0.,0.],[0.,0.]]
z=s.Symbol('z');poly=s.Poly(sum((-1)**r*s.binomial(11,r)*(center[r]+s.I*center[r+10])*z**r for r in range(10)),z)
points=[point(complex(r),1) for r in s.nroots(poly,n=30,maxsteps=1000)]+[point(None,2)]
states.append({'id':'gmax-n11-t3-certified','n':11,'achievedOrder':3,'label':'Certified algebraic state with one double point',
 'stateLatex':r'|\psi\rangle=\sum_{r=0}^{9}\sqrt{\binom{11}{r}}(x_r+i y_r)|D_{11}^{(r)}\rangle',
 'exactCoefficients':None,'coefficients':coeff,'support':list(range(10)),'partition':[2]+[1]*9,'maxMultiplicity':2,'distinctPoints':10,
 'majoranaPoints':points,'polynomialLatex':r'P(z)=\sum_{r=0}^{9}(-1)^r\binom{11}{r}(x_r+i y_r)z^r',
 'verification':{'method':'Exact rational Banach certificate of the unique root in a fixed-coordinate box; exact Sylvester-matrix bound proves nine simple finite roots and a double point at infinity.','exactAC':True,'exactMultiplicity':True,'certificatePath':'./research/n11-g2/witness.json'},
 'notes':'The coefficients are defined by the unique solution of sixteen rational quadratic equations in the specified box of radius 10⁻³⁰. Displayed decimals approximate that certified algebraic state. Four coordinates are fixed; uniqueness is within this slice, not among all states. The 2017 numerical vector is a starting point, not an exact witness by itself.',
 'references':[{'file':'PRA96_032304(2017).pdf','page':11,'location':'Appendix A, Eq. (A1), numerical seed; exact existence and multiplicities independently certified here.'}]})

proofs=[
 {'id':'moment-principle','title':'From a repeated point to exact moment bounds','summary':'A repeated Majorana point forces a truncated Dicke distribution in a suitable pole frame.','formulaLatex':r'K\in\{g,\ldots,n\},\qquad\mathbb E[K^a]=\frac1{n+1}\sum_{k=0}^{n}k^a\quad(0\le a\le t)',
  'bodyHtml':body(r'<p>Rotate a point of multiplicity \(g\) to the pole at which \(c_0=\cdots=c_{g-1}=0\). With \(\Pr(K=k)=|c_k|^2\), the diagonal anticoherence equations fix the first \(t\) moments to those of the uniform distribution on \(0,\ldots,n\).</p><p>If a polynomial \(p\) of degree at most \(t\) is nonnegative on every allowed integer \(g,\ldots,n\), but its uniform expectation is negative, that multiplicity is impossible. All displayed upper bounds have exact rational polynomial certificates. Feasible moment distributions alone do not establish anticoherence.</p><p><a href="./research/universal-bounds.json">Cell-by-cell upper bounds</a> · <a href="./research/moment-bounds.json">Full rational certificates</a> · <a href="./research/verify_moment_bounds.py">Verifier</a></p>'),
  'references':[{'file':'PRA92_052333(2015).pdf','page':4,'location':'Eqs. (16)–(17); exact moment-bound deductions supplied here.'}]},
 {'id':'order-one','title':'First order: the half-size bound is attained','summary':'The largest point multiplicity is the integer part of half the qubit number.','formulaLatex':r'g_{\max}(n,1)=\left\lfloor\frac n2\right\rfloor',
  'bodyHtml':body(r'<p>The mean \(\mathbb E[K]=n/2\) and support \(K\ge g\) imply \(g\le n/2\). Even sizes attain the bound with \(|D_n^{(n/2)}\rangle\), whose partition is \((n/2,n/2)\). Odd sizes attain it with \((|D_n^{(0)}\rangle+\sqrt n\,|D_n^{((n+1)/2)}\rangle)/\sqrt{n+1}\).</p><p>These are attaining representatives, not a classification of every state with the same maximal multiplicity.</p>'),
  'references':[{'file':'PRA96_032304(2017).pdf','page':9,'location':'Eqs. (58)–(59)'},{'file':'PRA90_032314(2014).pdf','page':4,'location':'Majorana-degeneracy nonexistence criterion'}]},
 {'id':'order-two','title':'Second order: a sharp one-third bound','summary':'The variance constraint gives a universal bound, attained by exact witnesses wherever second-order states exist in the core.','formulaLatex':r'g\le\left\lfloor\frac{n-1}{3}\right\rfloor',
  'bodyHtml':body(r'<p>On the allowed support, \((K-g)(n-K)\ge0\). Its expectation in a second-order anticoherent state is \(n(n-1-3g)/6\), giving the bound.</p><p>For every \(n\ge7\), put \(g=\lfloor(n-1)/3\rfloor\), \(\delta=n-3g-1\), and \(L=n-g\). The positive amplitudes on Dicke support \(0,3,L\) with probabilities \(p_0=[9(n-2g)-n\delta]/(18L)\), \(p_3=n\delta/[18(L-3)]\), and \(p_L=n[3(L-3)-\delta]/[6L(L-3)]\) attain the bound. The support gaps remove every off-diagonal condition through order two. When \(\delta=0\), the middle amplitude vanishes.</p><p>The four-qubit tetrahedral and six-qubit octahedral states cover the smaller admissible sizes. The sizes \(n=2,3,5\) admit no second-order state.</p><p><a href="./research/build_witnesses.py">Exact construction and reduced-matrix verifier</a></p>'),
  'references':[{'file':'PRA96_032304(2017).pdf','page':9,'location':'Eq. (60), the two-term subfamily; the all-size construction and sharp bound are independently derived here.'}]},
 {'id':'order-three','title':'Third order: a positive moment matrix','summary':'A two-by-two positive matrix limits how many points can coincide.','formulaLatex':r'g\le\left\lfloor\frac n2-\sqrt{\frac{n(n+2)}{12}}\right\rfloor',
  'bodyHtml':body(r'<p>Let \(X=K-n/2\), \(a=n/2-g\), and \(V=n(n+2)/12\). Third-order anticoherence gives \(\mathbb E[X]=\mathbb E[X^3]=0\) and \(\mathbb E[X^2]=V\). The nonnegative weight \(X+a\) makes the moment matrix \(\begin{pmatrix}a&V\\V&aV\end{pmatrix}\) positive semidefinite. Its determinant \(V(a^2-V)\) yields the bound.</p><p>At eleven qubits the attaining state is certified by a rational fixed-point argument. At twelve qubits an elementary three-term state attains two with partition \((2,2,2,2,2,2)\). Changing its middle phase gives the same maximal multiplicity with partition \((2,2,1,1,1,1,1,1,1,1)\).</p>'),
  'references':[{'file':'PRA92_052333(2015).pdf','page':4,'location':'Moment characterization; the matrix inequality is an independent catalogue deduction.'}]},
 {'id':'n12-higher','title':'Fourth and fifth order at twelve qubits: every point is simple','summary':'A short polynomial rules out even one double point.','formulaLatex':r'p(K)=(K-2)(12-K)(K-6)(K-7),\qquad\frac1{13}\sum_{k=0}^{12}p(k)=-42',
  'bodyHtml':body(r'<p>If a point had multiplicity at least two, the pole-frame support would be \(K\in\{2,\ldots,12\}\). The polynomial above is nonnegative at every allowed integer. Fourth-order anticoherence would instead force its expectation to be \(-42\), a contradiction. Thus \(g\le1\). Existing fourth- and fifth-order representatives attain one.</p>'),
  'references':[{'file':'PRA96_032304(2017).pdf','page':10,'location':'Eq. (62) supplies attaining states; the exclusion polynomial is independently derived.'}]},
 {'id':'n11-algebraic','title':'An exact eleven-qubit state with a double point','summary':'The numerical source vector is replaced by a uniquely specified algebraic state with a rigorous existence certificate.','formulaLatex':r'\rho_3=\frac{I_4}{4},\qquad\lambda=(2,1^9),\qquad g=2',
  'bodyHtml':body(r'<p>Write \(c_r=\sqrt{\binom{11}{r}}(x_r+i y_r)\) for \(0\le r\le9\), with \(c_{10}=c_{11}=0\). Four real coordinates are fixed rationally. The sixteen remaining real reduced-matrix equations have a unique root in the specified rational box of radius \(10^{-30}\).</p><p>An exact rational Banach bound proves that root exists. A separate inverse-perturbation bound for the polynomial’s Sylvester matrix proves that all nine finite roots are simple throughout the box. The point at infinity has multiplicity exactly two. Both arguments were independently checked.</p><p><a href="./research/n11-g2/proof.md">Full argument</a> · <a href="./research/n11-g2/formula.tex">LaTeX</a> · <a href="./research/n11-g2/witness.json">Algebraic definition and coefficient intervals</a> · <a href="./research/n11-g2/verify.py">Existence verifier</a> · <a href="./research/n11-g2/verify_majorana.py">Multiplicity verifier</a></p>'),
  'references':[{'file':'PRA96_032304(2017).pdf','page':11,'location':'Eq. (A1), numerical seed only. Exact certification is independently supplied here.'}]}
]
for p in ac['proofs']:
    q=dict(p);q['id']='ac-'+p['id'];q['bodyHtml']=q['bodyHtml'].replace('./research/reported-exclusions/','./research/ac-exclusions/');proofs.append(q)

# Additional states are defined by isolated algebraic parameters, with exact
# moment and full-multiplicity certificates in the accompanying constructor.
constructions=read('research/ac-order-constructions.json')
assert constructions['allVerifiedExactly']
construction_certificates={c['id']:c for c in constructions['certificates']}
for st in constructions['states']:
    if st['exactCoefficients'] is None:
        cert=construction_certificates[st['id']]
        lo,hi=map(s.Rational,cert['rootBracket'])
        definition={'polynomialLatex':cert['parameterPolynomialLatex']+'=0',
                    'intervalLatex':s.latex(lo)+'<k<'+s.latex(hi),
                    'recoveryLatex':'x='+cert['xRecoveryLatex']+r',\qquad u=\sqrt{x}>0'}
        states.append(dict(st,constructionProofId='partition-higher-order-constructions',algebraicDefinition=definition))
proofs.insert(0,{
 'id':'partition-higher-order-constructions',
 'title':'Nine improved lower bounds from exact second-order states',
 'summary':'Six previously unresolved partitions now have proven maximum two. Three eleven-qubit partitions improve from 1–3 to 2–3.',
 'formulaLatex':r'\sum_r(2r-n)|c_r|^2=0,\qquad\sum_r[3(2r-n)^2-n(n+2)]|c_r|^2=0',
 'bodyHtml':body(r'<p>Every new witness has Dicke support separated by at least three indices, so all off-diagonal entries of \(\rho_2\) vanish. The two exact moment identities above, together with normalization, give \(\rho_2=I_3/3\).</p>'
 r'<p>At eight qubits, probabilities \((7,4,16)/27\) on indices \((0,3,6)\), with phase \(i\) on the middle coefficient, give partition \((2,1^6)\). At ten qubits, positive amplitudes with probabilities \((1/6,1/3,1/2)\) on indices \((0,3,8)\) give \((2,1^8)\). Exact polynomial gcd calculations establish the remaining simple roots. Both partitions have proven maximum two.</p>'
 r'<p>The cyclic family \(P(z)=z^a(z^3-u)^2(z^3-ku)\), with \(u=\sqrt{x}>0\), adds \((2^3,1^3)\) at \(n=9\), \((2^3,1^4)\) at \(n=10\), \((2^3,1^5)\) and \((2^4,1^3)\) at \(n=11\), and \((3,2^3,1^3)\) at \(n=12\). A specified integer polynomial and rational isolating interval define \(k\); a rational function fixes \(x\). The certificate proves \(x>0\), \(k\notin\{0,1\}\), and zero remainder for both original moment equations. Thus the two cubics are squarefree and coprime; the poles give the remaining multiplicities. The family with cubic power three, \(P(z)=(z^3-u)^3(z^3-ku)\), also gives \((3^3,1^3)\) at \(n=12\), with its own isolated algebraic parameter. The nine-, ten- and twelve-qubit cases attain their rigorous upper bound of two.</p>'
 r'<p>A further eleven-qubit state uses probabilities \((61/126,1/36,41/84)\) on \((2,5,9)\). Its two double poles and seven simple finite roots give \((2^2,1^7)\). These three eleven-qubit partitions now satisfy \(2\le t_{\max}\le3\); no third-order exclusion is inferred.</p>'
 r'<p><a href="./research/ac-order-sparse.json">Sparse certificates</a> · <a href="./research/ac-order-constructions.json">Algebraic definitions and certificates</a> · <a href="./research/ac-order-improvements-audit.json">Independent audit</a> · <a href="./research/ac-order-constructions.py">Reproducible constructor</a></p>'),
 'references':[]})

proofs.append({
 'id':'n12-four-double-obstruction',
 'title':'Four double points at twelve qubits: a restricted exclusion',
 'summary':'Antipodal double-point pairs and rotational symmetry beyond a single half-turn axis are excluded at third order. The general partition remains open.',
 'formulaLatex':r'P=Q^2R,\quad Q=z^4+2hz^2+1,\quad -1<h<1\quad\Longrightarrow\quad t<3\ \text{for }\lambda=(2^4,1^4)',
 'bodyHtml':body(r'<p>A rotation puts any two distinct antipodal double-point axes into the stated quartic normal form. For arbitrary complex \(R\), third-order anticoherence forces \(R\) to be time-reversal real up to phase: the ten-dimensional odd part of its density matrix maps injectively into the three-qubit reduction. An exact determinant and real-root counts verify this throughout \(-1<h<1\).</p>'
 r'<p>The remaining real equations have only the projective solution \(P=z^2(z^4+1)^2\), whose full partition is \((2^6)\), including its double point at infinity. It lies outside the target partition. The proof therefore excludes paired double-point axes, not arbitrary doubled-point geometry.</p><p>Separate exact orbit and polynomial arguments also exclude threefold, fourfold and \(D_2\) symmetry. Any hypothetical third-order state in this partition must have proper rotational group \(C_1\) or \(C_2\). Those remaining cases have not been excluded.</p>'
 r'<p>The catalogue retains \(2\le t_{\max}(12,(2^4,1^4))\le3\). An explicit family with the required partition approaches a third-order state in \((2^6)\), so arbitrarily small numerical residuals alone cannot establish attainment.</p>'
 r'<p><a href="./research/n12-four-double/index.html">Complete proof, further restrictions and search results</a> · <a href="./research/n12-four-double/antipodal-obstruction-verification.json">Exact verification</a> · <a href="./research/n12-four-double/verify-antipodal-obstruction.py">Reproduce the proof</a></p>'),
 'references':[]})

bounds={(r['n'],r['t']):r for r in read('research/universal-bounds.json')['records']}
reports={r['n']:r['reportedGmax'] for r in read('research/findings.json')['figure6']['rows']}
entries=[]
for a in ac['entries']:
    n,t=a['n'],a['t']
    if n>12 or t>6:continue
    candidates=sorted([st for st in states if st['n']==n and st['achievedOrder']>=t],key=lambda st:(-st['maxMultiplicity'],-st['achievedOrder'],st['id']))
    entry={'id':f'n{n}-t{t}','n':n,'t':t,'sourceReport':reports[n].get(str(t)),'references':[{'file':'PRA96_032304(2017).pdf','page':9,'location':'Fig. 6, source numerical survey; maxima are independently certified here.'}],'stateIds':[st['id'] for st in candidates]}
    if a['status']=='impossible':
        entry.update(status='no_ac',value=None,lower=None,upper=None,proofId='ac-'+a['proofId'],boundLatex='',notes='No pure symmetric state has this anticoherence order at this qubit number. Maximal multiplicity is therefore undefined, rather than zero.')
        entry['references']=a['references']
    else:
        upper=bounds[n,t]['upperBoundOnMaxMultiplicity'];lower=max(st['maxMultiplicity'] for st in candidates)
        assert lower<=upper,(n,t,lower,upper)
        entry.update(status='proven' if lower==upper else 'bounded',value=lower if lower==upper else None,lower=lower,upper=upper,
                     proofId={1:'order-one',2:'order-two',3:'order-three',4:'n12-higher',5:'n12-higher'}[t],boundLatex=rf'g_{{\max}}({n},{t})\le {upper}',
                     notes='The maximum is attained by a verified state and matched by a rigorous upper bound.' if lower==upper else 'The catalogue supplies a verified lower bound and a rigorous upper bound.',
                     momentCertificate=bounds[n,t]['nextMultiplicityExclusion'])
    entries.append(entry)
sources=[source for source in ac['sources'] if source.get('id')!='mixed-ac-2026']
for source in sources:
    if source['filename']=='PRA96_032304(2017).pdf':source['notes']='Primary numerical degeneracy diagram (Fig. 6), exact families, and the numerical seed used for the certified eleven-qubit double-point state.'
    elif source['filename']=='PRA90_032314(2014).pdf':source['notes']='Majorana multiplicities, SLOCC degeneracy partitions, and the half-size nonexistence criterion.'
    elif source['filename']=='PRA92_052333(2015).pdf':source['notes']='Exact anticoherence equations, pole support restrictions, and the distinction between degeneracy partitions and SLOCC classes.'
    else:source['notes']='Supporting source retained from the companion AC catalogue; used for definitions, witnesses, or the independently audited existence table.'
counts=Counter(e['status'] for e in entries)
data={'schemaVersion':'1.0','title':'Majorana degeneracy of anticoherent states','compiledDate':'2026-10-28','numbers':list(range(2,13)),'orders':list(range(1,7)),
 'scope':'Pure symmetric states; largest point multiplicity among states of at least the selected anticoherence order. Complete n=2–12 core.',
 'entries':entries,'states':states,'proofs':proofs,'sources':sources,
 'conventions':{'polynomial':r'P(X,Y)=\sum_{r=0}^{n}(-1)^r\sqrt{\binom nr}\,c_r X^{n-r}Y^r','coordinates':'For z=Y/X a root corresponds to the spinor ratio b/a=1/z. Zero is the south pole and infinity is the north pole. Point locations are displayed numerically; multiplicities are exactly certified.','partition':'The full multiplicity partition sums to n. The maximal part is g; the number of parts is the number of distinct points. Neither maximal multiplicity nor the full partition determines a state or, generally, a SLOCC class.'},
 'validation':{'provenCells':counts['proven'],'noACCells':counts['no_ac'],'boundedCells':counts['bounded'],'stateCount':len(states),'allCoreCellsSettled':counts['bounded']==0,'radicalStateCount':sum(st['exactCoefficients'] is not None for st in states),'certifiedImplicitStateCount':sum(st['exactCoefficients'] is None for st in states)}}
assert len(entries)==56 and counts['proven']==26 and counts['no_ac']==30 and counts['bounded']==0
assert all(sum(st['partition'])==st['n'] and sum(p['multiplicity'] for p in st['majoranaPoints'])==st['n'] for st in states)
data=reindex(data)
data=attach_library(data,ROOT)
data=attach_acliege_exact(data,ROOT)
for name,value in [('catalogue.json',json.dumps(data,ensure_ascii=False,indent=2)+'\n'),('catalogue.js','window.MAJORANA_CATALOGUE='+json.dumps(data,ensure_ascii=False,separators=(',',':'))+';\n')]:
    (ROOT/'data'/name).write_text(value)
with (ROOT/'data/catalogue.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['n','partition','status','value','lower','upper','reportedOrder','coverage','attainment','proofId'],extrasaction='ignore');w.writeheader()
    w.writerows(dict(e,partition=','.join(map(str,e['partition']))) for e in data['entries'])
print('Built',len(data['entries']),'partition records,',len(states),'verified representatives:',data['validation'])
