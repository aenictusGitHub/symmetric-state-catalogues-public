#!/usr/bin/env python3
"""Exact adjacent-triple second-order AC construction.
All three bands have the SAME orientation. Reversing one band does not
permit triangle cancellation: the two directed J_+^2 entries are distinct.
"""
from pathlib import Path
import json
import sympy as s

ROOT = Path(__file__).resolve().parent

def construct(n, k, h, offset=0):
    assert n >= 4 and k >= 1 and h % 2 == 0 and h >= 2*k + 2
    rows = [[offset + 2*i + h*l for l in range(3)] for i in range(k)]
    assert rows[0][0] >= 0 and rows[-1][-1] <= n
    variance = s.Rational(n*(n+2), 12)
    weights = []
    for row in rows:
        m = [s.Rational(n, 2) - r for r in row]
        w = [s.factor((variance + m[(l+1)%3]*m[(l+2)%3]) /
                      ((m[l]-m[(l+1)%3])*(m[l]-m[(l+2)%3])))
             for l in range(3)]
        assert all(v > 0 for v in w)
        weights.append(w)
    phases = [[s.Integer(1)]*3]
    edges = []
    for i in range(k-1):
        A = [s.factor(weights[i][l]*weights[i+1][l]*(n-r)*(n-r-1)*(r+1)*(r+2))
             for l,r in enumerate(rows[i])]
        D = s.factor(2*(A[0]*A[1]+A[1]*A[2]+A[0]*A[2])-sum(x*x for x in A))
        assert D >= 0
        z1 = s.simplify((A[2]-A[0]-A[1]+s.I*s.sqrt(D))/(2*s.sqrt(A[0]*A[1])))
        z2 = s.simplify(-(s.sqrt(A[0])+s.sqrt(A[1])*z1)/s.sqrt(A[2]))
        z = [s.Integer(1), z1, z2]
        assert all(s.simplify(q*s.conjugate(q)-1)==0 for q in z)
        assert s.simplify(sum(s.sqrt(A[l])*z[l] for l in range(3))) == 0
        phases.append([s.simplify(phases[-1][l]*z[l]) for l in range(3)])
        edges.append({'squaredLengths':list(map(str,A)), 'discriminant':str(D),
                      'phaseIncrements':list(map(str,z))})
    B = s.zeros(n+1,k)
    for i in range(k):
        for l,r in enumerate(rows[i]):
            B[r,i] = s.simplify(s.sqrt(weights[i][l])*phases[i][l])
    return rows,weights,edges,B

def verify(n,B):
    k=B.cols;j=s.Rational(n,2)
    Jz=s.diag(*[j-r for r in range(n+1)])
    Jp=s.zeros(n+1)
    for r in range(1,n+1):Jp[r-1,r]=s.sqrt(r*(n-r+1))
    ops={'normalization':(s.eye(n+1),s.eye(k)),
         'Jz':(Jz,s.zeros(k)),
         'Jz_squared':(Jz**2,j*(j+1)/3*s.eye(k)),
         'Jplus':(Jp,s.zeros(k)),
         'Jplus_Jz':(Jp*Jz,s.zeros(k)),
         'Jplus_squared':(Jp**2,s.zeros(k))}
    checks={name:all(s.simplify(v)==0 for v in B.conjugate().T*op*B-target)
            for name,(op,target) in ops.items()}
    assert all(checks.values())
    return checks

if __name__=='__main__':
    n,k,h,o=22,3,8,1
    rows,weights,edges,B=construct(n,k,h,o)
    checks=verify(n,B)
    ident='adjacent-triples-n22-k3-t2'
    basis=[]
    for i in range(k):
        vec=list(B[:,i]);terms=[]
        for r in rows[i]:terms.append('('+s.latex(B[r,i])+')\\,|D_{22}^{('+str(r)+')}\\rangle')
        basis.append({'id':f'{ident}-basis-{i+1}','n':n,'label':f'Adjacent-triple exact basis vector {i+1}',
                      'support':rows[i], 'exactCoefficients':list(map(str,vec)),
                      'coefficients':[[float(s.re(v).evalf(17)),float(s.im(v).evalf(17))] for v in vec],
                      'stateLatex':f'|u_{i+1}\\rangle='+ '+'.join(terms)})
    compact=[r'|u_1\rangle=\frac12|D_{22}^{(1)}\rangle+\frac12|D_{22}^{(9)}\rangle+\frac1{\sqrt2}|D_{22}^{(17)}\rangle',
             r'|u_2\rangle=\sqrt{\frac{11}{32}}|D_{22}^{(3)}\rangle+\frac{\sqrt5}{4}z|D_{22}^{(11)}\rangle+\sqrt{\frac{11}{32}}w|D_{22}^{(19)}\rangle',
             r'|u_3\rangle=\frac1{\sqrt2}|D_{22}^{(5)}\rangle+\frac{w}{2}|D_{22}^{(13)}\rangle+\frac{w^2}{2}|D_{22}^{(21)}\rangle']
    definitions=r'z=\frac{-37+i\sqrt{9551}}{2\sqrt{2730}},\qquad w=\frac{-5-i\sqrt{9551}}{6\sqrt{266}},\qquad |z|=|w|=1.'
    for i,b in enumerate(basis):
        b['stateLatex']=compact[i]
        if i:b['coefficientDefinitionLatex']=definitions
    result={'id':ident,'n':n,'t':2,'dimension':k,'label':'Analytical adjacent-triple second-order AC subspace',
            'evidence':'exact','basis':basis,'parameters':{'spacing':h,'offset':o},
            'weights':[[str(x) for x in w] for w in weights],'edges':edges,
            'derivation':'Three equally oriented Dicke bands; the first two moments fix rational weights. Adjacent Jplus-squared cross terms cancel by exact triangle phases. All other cross terms vanish by support separation. This establishes K(22,2)>=3, not maximality.',
            'verification':checks}
    (ROOT/'adjacent-triples.json').write_text(json.dumps(result,indent=2)+'\n')
    (ROOT/'adjacent-triples-verification.json').write_text(json.dumps({'n':n,'dimension':k,'t':2,'exact':True,'checks':checks},indent=2)+'\n')
    print(json.dumps({'id':ident,'checks':checks}))
