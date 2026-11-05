#!/usr/bin/env python3
"""Construct and certify seven exact second-order AC witnesses.

Run with Python 3 and SymPy. The generated JSON contains exact defining
polynomials, rational isolating intervals, proof checks and display vectors.
Decimal coefficients are for visualization; the algebraic definitions are
what establish exact attainment and full Majorana multiplicity.
"""
import json
import math
from pathlib import Path
from fractions import Fraction
import sympy as s

HERE = Path(__file__).resolve().parent
OUT = HERE / 'ac-order-constructions.json'
x, k, z = s.symbols('x k z', real=True)


def integer_polynomial(expr, *variables):
    return s.Poly(expr, *variables).clear_denoms()[1].primitive()[1]


def interval_multiply(a, b):
    products = [u * v for u in a for v in b]
    return min(products), max(products)


def interval_evaluate(poly, interval):
    out = (s.Rational(0), s.Rational(0))
    for coefficient in s.Poly(poly, k).all_coeffs():
        out = interval_multiply(out, interval)
        out = out[0] + coefficient, out[1] + coefficient
    return out


def root_point(value, multiplicity):
    value = complex(value)
    denominator = 1 + abs(value) ** 2
    return {'x': 2 * value.real / denominator,
            'y': -2 * value.imag / denominator,
            'z': (abs(value) ** 2 - 1) / denominator,
            'multiplicity': multiplicity}


def cubic_roots(value):
    value = complex(value)
    radius, angle = abs(value) ** (1 / 3), math.atan2(value.imag, value.real)
    return [radius * complex(math.cos((angle + 2 * math.pi * j) / 3),
                             math.sin((angle + 2 * math.pi * j) / 3))
            for j in range(3)]


def make_algebraic(n, a, b, root_bracket, repeated_multiplicity=2):
    m = repeated_multiplicity
    assert n == a + 3*(m+1) + b
    support = [a + 3 * i for i in range(m+2)]
    # P=z^a(z^3-u)^m(z^3-ku), with x=u^2>0.
    coefficient_factors = [s.binomial(m,i-1)+k*s.binomial(m,i) for i in range(m+2)]
    squared_polynomial_coefficients = [v*v*x**(m+1-i) for i,v in enumerate(coefficient_factors)]
    weights = [v/s.binomial(n, r) for r, v in zip(support, squared_polynomial_coefficients)]
    F = integer_polynomial(sum((2*r-n)*w for r,w in zip(support,weights)), x, k).as_expr()
    G = integer_polynomial(sum((3*(2*r-n)**2-n*(n+2))*w for r,w in zip(support,weights)), x, k).as_expr()
    resultant = s.Poly(s.resultant(F, G, x), k)
    lo, hi = map(s.Rational, root_bracket)
    factors = s.factor_list(resultant)[1]
    matches = [factor for factor, _ in factors if factor.count_roots(lo, hi) == 1]
    assert len(matches) == 1
    R = matches[0].monic().clear_denoms()[1].primitive()[1]
    if R.LC() < 0: R = -R
    assert R.is_irreducible and R.count_roots(lo, hi) == 1
    isolated = R.refine_root(lo, hi, eps=s.Rational(1, 10**35))
    assert R.count_roots(*isolated) == 1
    # A linear subresultant fixes x rationally; all checks are modulo R.
    linear = next(q for q in reversed(s.subresultants(F, G, x)) if s.degree(q,x) == 1)
    recovery = s.cancel(-linear.subs(x, 0) / s.diff(linear, x))
    numerator, denominator = map(s.expand, s.fraction(recovery))
    numerator = s.rem(s.Poly(numerator,k),R).as_expr()
    denominator = s.rem(s.Poly(denominator,k),R).as_expr()
    recovery = s.cancel(numerator / denominator)
    numerator, denominator = map(s.expand,s.fraction(recovery))
    assert s.gcd(s.Poly(denominator,k),R).degree() == 0
    remainders = []
    for equation in [F,G]:
        exact_numerator = s.fraction(s.cancel(equation.subs(x,recovery)))[0]
        rem = s.rem(s.Poly(exact_numerator,k),R)
        assert rem.is_zero
        remainders.append(str(rem.as_expr()))
    ni = interval_evaluate(numerator, isolated)
    di = interval_evaluate(denominator, isolated)
    assert (ni[0]>0 and di[0]>0) or (ni[1]<0 and di[1]<0)
    xi = interval_multiply(ni,(1/di[1],1/di[0]))
    assert xi[0] > 0
    assert not (isolated[0] <= 0 <= isolated[1])
    assert not (isolated[0] <= 1 <= isolated[1])
    # High precision is used only for the display record. Isolation above is exact.
    kv = s.nsolve(R.as_expr(), k, (isolated[0]+isolated[1])/2,
                  prec=90, maxsteps=100, verify=True)
    assert isolated[0] < kv < isolated[1]
    xv = s.N(recovery.subs(k,kv),85)
    uv = s.sqrt(xv)
    raw = [(-1)**(m+1-i)*v.subs(k,kv)*uv**(m+1-i) for i,v in enumerate(coefficient_factors)]
    vector = [s.Float(0,85)]*(n+1)
    for r, v in zip(support,raw): vector[r] = (-1)**r*v/s.sqrt(s.binomial(n,r))
    norm = s.sqrt(sum(v*v for v in vector))
    vector = [s.N(v/norm,75) for v in vector]
    points = ([{'x':0.0,'y':0.0,'z':-1.0,'multiplicity':a}] if a else [])
    points += [root_point(q,m) for q in cubic_roots(float(uv))]
    points += [root_point(q,1) for q in cubic_roots(float(kv*uv))]
    if b: points += [{'x':0.0,'y':0.0,'z':1.0,'multiplicity':b}]
    partition = sorted(([a] if a else [])+[m]*3+[1]*3+([b] if b else []),reverse=True)
    ident = f'ac-order-c3-n{n}-'+(f'm{m}-' if m!=2 else '')+f'a{a}-b{b}'
    bracket_text = f'{root_bracket[0]}<k<{root_bracket[1]}'
    polynomial = ('z^{'+str(a)+'}' if a else '')+r'(z^3-u)^{'+str(m)+r'}(z^3-ku)'
    state = {
        'id':ident,'n':n,'achievedOrder':2,'witnessesOrders':[1,2],
        'label':'Exact cyclic algebraic construction',
        'status':'proven-existence','support':support,'partition':partition,
        'maxMultiplicity':max(partition),'distinctPoints':len(partition),
        'stateLatex':r'|\psi\rangle=\frac{1}{\sqrt{Z}}\sum_{r\in\{'+','.join(map(str,support))+r'\}}\frac{(-1)^r[z^r]P(z)}{\sqrt{\binom{'+str(n)+r'}{r}}}|D_{'+str(n)+r'}^{(r)}\rangle,\quad Z=\sum_r\frac{|[z^r]P(z)|^2}{\binom{'+str(n)+r'}{r}}',
        'polynomialLatex':r'P(z)='+polynomial+r',\quad u=\sqrt{x}>0',
        'exactCoefficients':None,
        'coefficients':[[float(v),0.0] for v in vector],
        'highPrecisionCoefficients':[[str(s.N(v,60)),'0'] for v in vector],
        'majoranaPoints':points,
        'normalization':'Normalized exactly by the positive sum Z; decimal display rounded.',
        'coefficientConvention':'P(z)=sum_r (-1)^r sqrt(binomial(n,r)) c_r z^r, up to its nonzero normalization.',
        'notes':f'An exact algebraic state: {bracket_text}; k is the unique root of the defining polynomial in the certificate, x is its stated rational function, and u=sqrt(x)>0. The algebraic definition, rather than rounded display coefficients, supplies the exact witness.',
        'verification':{'method':'Exact rational root isolation, positive rational recovery, vanishing moment-polynomial remainders, and explicit coprime squarefree cubic factors.',
                        'exactAC':True,'exactMultiplicity':True,
                        'certificatePath':'./research/ac-order-constructions.json'},
        'references':[{'title':'Exact cyclic construction certificate','url':'./research/ac-order-constructions.json'}],
        'provenance':'New exact algebraic construction verified in this catalogue.'}
    cert = {
        'id':ident,'n':n,'poleAtZero':a,'poleAtInfinity':b,'repeatedMultiplicity':m,
        'support':support,'targetOrder':2,'partition':partition,
        'parameterPolynomial':str(R.as_expr()),
        'parameterPolynomialLatex':s.latex(R.as_expr()),
        'rootBracket':[str(lo),str(hi)],'rootInterval':list(map(str,isolated)),
        'rootCountInInterval':1,'irreducible':True,
        'xNumerator':str(numerator),'xDenominator':str(denominator),
        'xPositiveInterval':list(map(str,xi)),
        'xRecoveryLatex':s.latex(recovery),
        'momentPolynomials':[str(F),str(G)],'exactRemainders':remainders,
        'parametersDecimal':{'k':str(s.N(kv,70)),'x':str(s.N(xv,70)),'u':str(s.N(uv,70))},
        'normalizationPositive':True,'offDiagonalRho2VanishesBySupportGaps':True,
        'squarefreeCubicFactors':'u is nonzero; z^3-u and z^3-ku have three simple nonzero roots each.',
        'factorsCoprime':'k is neither zero nor one, so the two cubic factors and the pole factor have disjoint roots.',
        'extraneousResultantRootsExcluded':'The chosen root is isolated, the x recovery denominator is coprime to its polynomial, x is positive, and both unsquared moment equations vanish modulo its polynomial. No squared elimination condition alone is used.'}
    return state, cert


def make_sparse():
    n=11
    probs={2:s.Rational(61,126),5:s.Rational(1,36),9:s.Rational(41,84)}
    assert sum(probs.values())==1
    assert sum((2*r-n)*p for r,p in probs.items())==0
    assert sum((3*(2*r-n)**2-n*(n+2))*p for r,p in probs.items())==0
    coeffs=[s.sqrt(probs.get(r,0)) for r in range(n+1)]
    P=s.Poly(sum((-1)**r*s.sqrt(s.binomial(n,r))*c*z**r for r,c in enumerate(coeffs)),z)
    quotient=s.expand(P.as_expr()/z**2)
    A,B,C=[quotient.coeff(z,i) for i in [0,3,7]]
    # A+B z^3+C z^7 has no repeated complex root: a common root would
    # give z^12=(-3B/7C)^3=(-7A/4B)^4, but these real quantities
    # have opposite signs. It has no zero root since A>0.
    lhs=s.simplify((-3*B/(7*C))**3);rhs=s.simplify((-7*A/(4*B))**4)
    assert A>0 and B<0 and C<0 and lhs<0 and rhs>0
    roots=s.nroots(quotient,n=60,maxsteps=200)
    points=[{'x':0.0,'y':0.0,'z':-1.0,'multiplicity':2}]
    points += [root_point(complex(r),1) for r in roots]
    points += [{'x':0.0,'y':0.0,'z':1.0,'multiplicity':2}]
    partition=[2,2]+[1]*7
    ident='ac-order-sparse-n11-2x2-1x7'
    state={'id':ident,'n':11,'achievedOrder':2,'witnessesOrders':[1,2],
           'label':'Exact sparse two-pole construction','status':'proven-existence',
           'support':list(probs),'partition':partition,'maxMultiplicity':2,'distinctPoints':9,
           'stateLatex':r'\sqrt{\frac{61}{126}}|D_{11}^{(2)}\rangle+\frac16|D_{11}^{(5)}\rangle+\sqrt{\frac{41}{84}}|D_{11}^{(9)}\rangle',
           'polynomialLatex':s.latex(P.as_expr()),
           'exactCoefficients':list(map(str,coeffs)),
           'coefficients':[[float(c),0.0] for c in coeffs],
           'highPrecisionCoefficients':[[str(s.N(c,60)),'0'] for c in coeffs],
           'majoranaPoints':points,'normalization':'sum_r |c_r|^2=1, checked exactly',
           'verification':{'method':'Exact rational moments and support gaps; exact sign contradiction excludes repeated roots of the remaining degree-seven factor.','exactAC':True,'exactMultiplicity':True,'certificatePath':'./research/ac-order-constructions.json'},
           'references':[{'title':'Exact sparse construction certificate','url':'./research/ac-order-constructions.json'}],
           'notes':'Two double poles and seven simple finite roots. This witnesses order two; the partition maximum may be higher.',
           'provenance':'New exact sparse construction verified in this catalogue.'}
    certificate={'id':ident,'n':11,'partition':partition,'probabilities':{str(r):str(p) for r,p in probs.items()},
                 'normalizationVerifiedExactly':True,'firstMomentVerifiedExactly':True,'secondMomentVerifiedExactly':True,
                 'offDiagonalRho2VanishesBySupportGaps':True,
                 'remainingFactor':str(quotient),'repeatedRootContradiction':{'lhs':str(lhs),'rhs':str(rhs),'lhsNegative':True,'rhsPositive':True}}
    return state,certificate


def main():
    cases=[(9,0,0,('72.68','72.69')),
           (10,0,1,('24.63','24.64')),
           (11,1,1,('10.11','10.12')),
           (11,0,2,('9.49','9.50')),
           (12,0,3,('-7.16','-7.15')),
           (12,0,0,('-6760.50','-6760.49'),3)]
    states=[];certificates=[]
    for case in cases:
        state,certificate=make_algebraic(*case)
        states.append(state);certificates.append(certificate)
    state,certificate=make_sparse();states.append(state);certificates.append(certificate)
    proof={
      'id':'partition-cyclic-algebraic-order2',
      'title':'Exact second-order constructions with repeated Majorana roots',
      'summary':'Six cyclic algebraic states and one sparse state establish second-order anticoherence while preserving their stated full multiplicity partitions.',
      'formulaLatex':r'\sum_r(2r-n)|c_r|^2=0,\qquad\sum_r\left[3(2r-n)^2-n(n+2)\right]|c_r|^2=0',
      'bodyHtml':'<p>All nonzero Dicke coefficients are separated by at least three indices. Consequently the off-diagonal entries of the two-qubit reduced state vanish. The two exact moment identities above then make every diagonal entry equal to one third.</p><p>The cyclic family uses P(z)=z^a(z^3−u)^m(z^3−ku), with m equal to two or three and u positive. A univariate integer polynomial and a rational isolating interval define k; a rational function fixes x=u². The constructor verifies that x is strictly positive and that both original moment equations have zero remainder modulo the defining polynomial. Since k is neither zero nor one, the cubic factors have distinct nonzero roots and are coprime. Their multiplicities, together with the projective poles, therefore give the full partition exactly.</p><p>The additional eleven-qubit sparse state has two double poles. Its remaining factor A+Bz³+Cz⁷ is squarefree: a repeated root would require (−3B/7C)³=(−7A/4B)⁴, whose two sides have opposite signs. Exact existence does not by itself prove that order two is the largest possible order for a partition.</p>',
      'references':[{'title':'Reproducible constructor and exact certificate','url':'./research/ac-order-constructions.json'}]}
    output={'schemaVersion':'1.0','allVerifiedExactly':True,'stateCount':len(states),
            'method':'Exact algebraic moment identities and full root multiplicities; decimals are display-only.',
            'states':states,'certificates':certificates,'proof':proof}
    OUT.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'allVerifiedExactly':True,'stateCount':len(states),'output':str(OUT),
                      'partitions':[{'n':v['n'],'partition':v['partition'],'order':2} for v in states]},indent=2))


if __name__=='__main__':main()
