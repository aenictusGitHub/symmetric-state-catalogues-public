#!/usr/bin/env python3
"""Rigorous inertia obstructions inside the specific catalogue t=1 parents.

All arithmetic used to certify pivot signs is rational interval arithmetic.
The result is not an upper bound on arbitrary second-order AC subspaces.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import json
import sympy as s
ROOT = Path(__file__).resolve().parent
SCALE = 10**45

def floor_fraction(x): return x.numerator // x.denominator
def ceil_fraction(x): return -floor_fraction(-x)
class I:
    def __init__(self,lo,hi=None):
        self.lo=F(lo); self.hi=F(lo if hi is None else hi)
    @classmethod
    def outward(cls,lo,hi):
        return cls(F(floor_fraction(lo*SCALE),SCALE), F(ceil_fraction(hi*SCALE),SCALE))
    def __add__(a,b):
        if not isinstance(b,I): b=I(b)
        return I.outward(a.lo+b.lo,a.hi+b.hi)
    __radd__=__add__
    def __neg__(a): return I(-a.hi,-a.lo)
    def __sub__(a,b): return a+-as_i(b)
    def __rsub__(a,b): return as_i(b)+-a
    def __mul__(a,b):
        b=as_i(b); p=[a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi]
        return I.outward(min(p),max(p))
    __rmul__=__mul__
    def __truediv__(a,b):
        b=as_i(b); assert b.lo*b.hi>0, 'A divisor interval contains zero'
        return a*I.outward(1/b.hi,1/b.lo)
    def __rtruediv__(a,b): return as_i(b)/a

def as_i(x): return x if isinstance(x,I) else I(x)
def evaluate(expr):
    if expr.is_Rational: return I(F(int(expr.p),int(expr.q)))
    if expr.is_Add: return sum((evaluate(a) for a in expr.args),I(0))
    if expr.is_Mul:
        result=I(1)
        for a in expr.args: result=result*evaluate(a)
        return result
    if expr.is_Pow:
        b,e=expr.args
        if e==s.Rational(1,2):
            assert b.is_Rational and b>=0
            q=F(int(b.p),int(b.q)); n=isqrt(floor_fraction(q*SCALE*SCALE))
            return I(F(n,SCALE),F(n+1,SCALE))
        if e==-1: return I(1)/evaluate(b)
        if e.is_Integer and e>=0:
            result=I(1)
            for _ in range(int(e)): result=result*evaluate(b)
            return result
    raise ValueError(str(expr))

def check(n,coefficients):
    rows=json.loads((ROOT.parent/'analytical-constructions/t1-bases.json').read_text())['subspaces']
    row=next(x for x in rows if x['n']==n and x['construction']['type']=='reflection-pairs')
    U=s.Matrix.hstack(*[s.Matrix([s.sympify(v) for v in b['exactCoefficients']]) for b in row['basis']])
    j=s.Rational(n,2); c=j*(j+1)/3
    Z=s.diag(*[j-r for r in range(n+1)]); raising=s.zeros(n+1)
    for r in range(1,n+1): raising[r-1,r]=s.sqrt(r*(n-r+1))
    X=(raising+raising.T)/2; Y=(raising-raising.T)/(2*s.I)
    assert (U.T*U).applyfunc(s.simplify)==s.eye(U.cols)
    for J in [X,Y,Z]: assert (U.H*J*U).applyfunc(s.simplify)==s.zeros(U.cols)
    A=coefficients[0]*(Z*Z-c*s.eye(n+1))+coefficients[1]*(X*X-Y*Y)+coefficients[2]*(X*Z+Z*X)
    M=(U.H*A*U).applyfunc(s.simplify)
    assert M==M.T
    intervals=[[evaluate(M[i,q]) for q in range(M.cols)] for i in range(M.rows)]
    L=[[I(0) for _ in range(M.cols)] for _ in range(M.rows)];D=[]
    for a in range(M.rows):
        L[a][a]=I(1)
        D.append(intervals[a][a]-sum((L[a][b]*L[a][b]*D[b] for b in range(a)),I(0)))
        assert D[-1].lo*D[-1].hi>0, 'Inconclusive pivot'
        for r in range(a+1,M.rows):
            L[r][a]=(intervals[r][a]-sum((L[r][b]*L[a][b]*D[b] for b in range(a)),I(0)))/D[a]
    signs=[1 if d.lo>0 else -1 for d in D]
    integer_boxes={18:[(-31,-29),(15,16),(-4,-3),(91,92),(54,55),(40,41),(5,6)],
      24:[(-53,-51),(231,232),(-34,-33),(107,108),(52,53),(-33,-32),(377,378),(78,79),(26,27)]}[n]
    assert all(a<d.lo<=d.hi<b for (a,b),d in zip(integer_boxes,D))
    return {'N':n,'parentId':row['id'],'parentDimension':U.cols,
      'quadrupoleCoefficients':dict(zip(['Jz^2-cI','Jx^2-Jy^2','{Jx,Jz}'],coefficients)),
      'pivotSigns':signs,'simpleIntegerPivotEnclosures':integer_boxes,'pivotRationalIntervals':[[str(d.lo),str(d.hi)] for d in D],
      'inertia':{'positive':signs.count(1),'negative':signs.count(-1),'zero':0},
      'nestedT2DimensionUpperBound':min(signs.count(1),signs.count(-1)),
      'scope':'Only subspaces of this explicit t=1 parent (and common collective rotations of it)',
      'verified':True}

if __name__=='__main__':
    data={'arithmetic':'Rational intervals with outward rounding to 45 decimal places at every operation',
      'proof':'Every t=2 child has zero compression of each traceless quadrupole. Sylvester inertia bounds its dimension by min(n_positive,n_negative)+n_zero.',
      'records':[check(18,[1,1,1]),check(24,[1,2,1])],
      'unrestrictedUpperBound':False}
    # 29/(3sqrt(3))>83/15 follows from the following positive squared gap.
    assert 145**2-3*83**2==358
    finite=[]
    for n in range(22,29):
        j=s.Rational(n,2); r=(s.sqrt(1+4*j*(j+1)/3)-1)/2
        bound=int(2*s.floor((j+2-r)/3)); target=min((n-2)//4,(n+5)//5)
        assert bound<target
        finite.append({'N':n,'universalParentBound':bound,'generalUpperBound':target})
    data['generalParentObstruction']={'minimumN':22,'asymptoticBound':'(1-1/sqrt(3))*N/3 + O(1)',
      'exactBoundForNge15':'2 floor((j+2-(sqrt(1+4j(j+1)/3)-1)/2)/3)',
      'squaredGapForAllNge29':358,'finiteN22to28':finite,
      'scope':'Only the explicitly constructed maximal first-order parents, not all first-order AC subspaces'}
    (ROOT/'reflection-inertia-verification.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({str(r['N']):{'inertia':r['inertia'],'bound':r['nestedT2DimensionUpperBound']}for r in data['records']},indent=2))
