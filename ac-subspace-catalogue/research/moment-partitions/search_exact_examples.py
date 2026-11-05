#!/usr/bin/env python3
"""Exact simple AC2 examples from moment partitions and joint compression.
Reconstructs all ordered two-qubit cross reductions; no numerical tolerance.
"""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parent

def cross_reduction(n,u,v):
    result=s.zeros(3)
    for a in range(3):
        for b in range(3):
            for ell in range(n-1):
                x,y=u[ell+a],v[ell+b]
                if x and y:
                    factor=s.sqrt(s.binomial(2,a)*s.binomial(2,b))*s.binomial(n-2,ell)/s.sqrt(s.binomial(n,ell+a)*s.binomial(n,ell+b))
                    result[a,b]+=x*s.conjugate(y)*factor
    return result.applyfunc(s.simplify)

def verify(n,basis):
    k=len(basis); checks=[]
    for i,u in enumerate(basis):
        for j,v in enumerate(basis):
            target=s.eye(3)/3 if i==j else s.zeros(3)
            assert cross_reduction(n,u,v)==target,(n,i,j)
            checks.append({'i':i,'j':j,'allEntriesExactlyZero':True})
    return {'orthonormal':True,'allCrossReductionsVerified':True,'exactArithmetic':True,'checks':checks,'scalarEntriesChecked':9*k*k}

def build_record(n,ident,label,vectors,derivation):
    basis=[]
    for a,vec in enumerate(vectors):
        support=[r for r,x in enumerate(vec) if x]
        terms=[s.latex(vec[r])+r'\,|D_{'+str(n)+r'}^{('+str(r)+r')}\rangle' for r in support]
        basis.append({'id':f'{ident}-basis-{a+1}','n':n,'label':f'Exact moment-construction basis vector {a+1}', 'support':support,
          'exactCoefficients':[str(s.simplify(x)) for x in vec],
          'coefficients':[[float(s.re(x).evalf(17)),float(s.im(x).evalf(17))] for x in vec],
          'stateLatex':r'|u_{'+str(a+1)+r'}\rangle='+'+'.join(terms)})
    return {'id':ident,'n':n,'t':2,'order':2,'dimension':len(vectors),'label':label,'evidence':'exact','basis':basis,'verification':verify(n,vectors),'derivation':derivation}

def main():
    a,b,v=s.symbols('a b v',positive=True)
    p=(v-b*b)/(a*a-b*b);q=(a*a-v)/(a*a-b*b)
    A=v*q/(a*a);B=v*p/(b*b);C=1-v*v/(a*a*b*b)
    symbolic={
       'oddNormalization':s.factor(p+q-1),
       'oddSecondMoment':s.factor(a*a*p+b*b*q-v),
       'evenNormalization':s.factor(A+B+C-1),
       'evenSecondMoment':s.factor(a*a*A+b*b*B-v),
       'mixedFirstMomentSquaredCancellation':s.factor(a*a*A*p-b*b*B*q)}
    assert all(x==0 for x in symbolic.values())
    u=s.zeros(13,1);w=s.zeros(13,1)
    u[0]=u[12]=s.sqrt(s.Rational(77,486));u[3]=u[9]=s.sqrt(s.Rational(35,243));u[6]=4*s.sqrt(2)/9
    w[0]=s.sqrt(s.Rational(5,54));w[12]=-w[0];w[3]=-s.sqrt(s.Rational(11,27));w[9]=-w[3]
    r12=build_record(12,'moment-compression-n12-k2-t2','Real radical basis from a five-level moment compression',[list(u),list(w)],
       'On m=+/-6,+/-3,0, combine even and odd reflection vectors. Their second moment is 14 and their mixed first moment cancels exactly. The global Dicke spacing is three, so all rank-one and rank-two shift terms vanish. The supports overlap: this is a joint compression, not a disjoint convex partition. It attains the already-known dimension two.')
    r12['basis'][0]['stateLatex']=r'\begin{aligned}|u_1\rangle&=\sqrt{\frac{77}{486}}\bigl(|D_{12}^{(0)}\rangle+|D_{12}^{(12)}\rangle\bigr)\\&\quad+\sqrt{\frac{35}{243}}\bigl(|D_{12}^{(3)}\rangle+|D_{12}^{(9)}\rangle\bigr)+\frac{4\sqrt2}{9}|D_{12}^{(6)}\rangle.\end{aligned}'
    r12['basis'][1]['stateLatex']=r'\begin{aligned}|u_2\rangle&=\sqrt{\frac5{54}}\bigl(|D_{12}^{(0)}\rangle-|D_{12}^{(12)}\rangle\bigr)\\&\quad-\sqrt{\frac{11}{27}}\bigl(|D_{12}^{(3)}\rangle-|D_{12}^{(9)}\rangle\bigr).\end{aligned}'
    vectors=[]
    for support,weights in [([0,15],[s.Rational(4,15),s.Rational(11,15)]),([7,22],[s.Rational(11,15),s.Rational(4,15)]),([3,10,18],[s.Rational(37,105),s.Rational(3,14),s.Rational(13,30)])]:
        vec=[s.Integer(0)]*23
        for r,p in zip(support,weights):vec[r]=s.sqrt(p)
        vectors.append(vec)
    r22=build_record(22,'moment-partition-n22-k3-t2','Real radical basis from two moment pairs and one triangle',vectors,
       'Partition the globally three-separated Dicke support {0,3,7,10,15,18,22} into {0,15}, {7,22}, {3,10,18}. The displayed positive rational probabilities have first spin moment zero and second spin moment 44. Orthogonality and all non-diagonal rank-one and rank-two conditions follow from the disjoint supports and global spacing. This simplifies the earlier dimension-three witness and does not increase the catalogue lower bound.')
    out={'schemaVersion':1,'scope':'Exact analytical witnesses; no new unrestricted maximum or improved lower bound is asserted.','subspaces':[r12,r22]}
    (ROOT/'search-exact-examples.json').write_text(json.dumps(out,indent=2)+'\n')
    report={'generalFiveLevelFormula':{'domain':'a>b>0, b^2<=v<=a*b; symmetric nodes +/-a,+/-b,0','p':str(p),'q':str(q),'A':str(A),'B':str(B),'C':str(C),'symbolicIdentities':{key:str(value) for key,value in symbolic.items()},'positivity':'Follows directly from b^2<=v<=a*b<a^2.','phaseSigns':'The odd b-pair has a minus sign; all other displayed amplitudes are positive.'},'examples':[{'n':r['n'],'k':r['dimension'],**r['verification']} for r in [r12,r22]],'passed':True}
    (ROOT/'search-exact-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':True,'examples':[(r['n'],r['dimension']) for r in [r12,r22]],'orderedCrossEntries':sum(r['verification']['scalarEntriesChecked'] for r in [r12,r22])}))
if __name__=='__main__':main()
