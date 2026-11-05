"""Combine the radical N24 witness and certified algebraic N26 record.

Requires mpmath to render the rational center as high-precision decimal values.
Those decimal values are never used as proof premises.
"""
from pathlib import Path
from fractions import Fraction
import copy,json
import mpmath as mp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def main():
    n24=json.loads((HERE/'n24-t7/states.json').read_text())['states'][0]
    n24=copy.deepcopy(n24)
    n24['stateLatex']=r'\begin{gathered}'+n24['stateLatex']+r'\\[6pt]'+n24['parametersLatex']+r'\end{gathered}'
    n24['references']=[{'title':'Twenty-four qubits: exact maximum order seven','url':'./research/n24-n26/index.html#n24'},
        {'title':'Analytical construction and universal upper-bound proof','url':'./research/n24-n26/n24-t7/proof.tex'},
        {'title':'Independent exact verification','url':'./research/n24-n26/n24-t7/independent-verification.json'}]
    n24['verification']['globalMaximumOrderProven']=True
    n24['verification']['maximumOrderUpperBound']=7
    n24['verification']['reductions']=[{'t':c['order'],'identityVerified':c['equalsSymmetricIdentity'],
        'frobeniusDistanceSquared':c['squaredFrobeniusDefect'],
        'purity':str(Fraction(c['squaredFrobeniusDefect'])+Fraction(1,c['order']+1))}
        for c in n24['verification']['reductions']]
    n24['verification']['nextOrderWitness']['differenceLatex']=r'\frac{13}{155040}'
    n24['majoranaCoordinatePolicy']='numeric-roots-retain-explicit-poles'
    n24['notes']='Exactly seventh-order AC. The rectangular spin-tensor argument excludes order eight for every pure symmetric N=24 state, proving the global maximum seven.'

    source=json.loads((ROOT/'data/catalogue.json').read_text())
    original=next(s for s in source['states'] if s['id']=='liege-acs.26.6')
    cert=json.loads((HERE/'n26-t6/certificate.json').read_text())
    center=list(map(Fraction,cert['center']))
    scales=cert['scaleFactors'];mp.mp.dps=150
    approximations=[]
    for r in range(27):
        factor=mp.sqrt(mp.binomial(26,r))/scales[r]
        approximations.append([mp.nstr(factor*(mp.mpf(center[k].numerator)/center[k].denominator),100)
            for k in [r,27+r]])
    n26={'id':'certified-algebraic-n26-t6','n':26,'achievedOrder':6,
        'witnessesOrders':list(range(1,7)),'label':'Certified algebraic state · order 6',
        'status':'proven-existence','coefficientEvidence':'implicit-algebraic',
        'stateLatex':r'|\psi\rangle=\sum_{r=0}^{26}\frac{\sqrt{\binom{26}{r}}}{d_r}\bigl(x_r^*+ix_{27+r}^*\bigr)|D_{26}^{(r)}\rangle',
        'support':list(range(27)),
        'coefficients':[[float(x),float(y)] for x,y in approximations],
        'highPrecisionCoefficients':approximations,
        'majoranaPoints':copy.deepcopy(original['majoranaPoints']),
        'partition':copy.deepcopy(original['partition']),
        'majoranaCoordinatePolicy':'Numeric rendering copied from the nearby refined Liège state; the exact state is defined by the certified polynomial root.',
        'normalization':'Exact rho_6 = I_7/7 implies sum_r |c_r|^2 = 1.',
        'coefficientConvention':'Normalized Dicke basis; r counts excitations. Floating components render the rational certificate center.',
        'implicitAlgebraicDefinition':{
            'certificateUrl':'./research/n24-n26/n26-t6/certificate.json',
            'proofUrl':'./research/n24-n26/n26-t6/proof.tex',
            'equations':'49 real rational quadratics equivalent to rho_6=I_7/7 after the stated binomial scaling.',
            'variableConvention':cert['variableConvention'],
            'activeVariables':cert['activeVariables'],'fixedVariables':cert['fixedVariables'],
            'radius':cert['radius'],
            'rootDefinition':'Unique real zero in the active-coordinate infinity-norm box around the saved rational center, with the five fixed coordinates equal to their saved rational values.',
            'centerCoefficientErrorBound':'2e-65 per real or imaginary Dicke component, relative to the exact rational-center conversion; ordinary floating values have additional rounding.',
            'closedFormRadicalsClaimed':False},
        'references':[{'title':'Twenty-six qubits: certified maximum order six','url':'./research/n24-n26/index.html#n26'},
            {'title':'Exact algebraic existence proof','url':'./research/n24-n26/n26-t6/proof.tex'},
            {'title':'Rational root certificate','url':'./research/n24-n26/n26-t6/certificate.json'},
            {'title':'Independent rational verification','url':'./research/n24-n26/n26-t6/independent-verification.json'},
            {'title':'Universal exclusion of seventh order','url':'./research/n24-n26/n26-t7/proof.tex'}],
        'provenance':'Certified exact algebraic zero near the refined Liège acs.26.6 numerical state; exact rational Banach contraction certificate.',
        'derivedFromStateId':'liege-acs.26.6','sourceLibraryRecordId':'acs.26.6',
        'notes':'Exact sixth-order existence is proved by an isolating rational contraction certificate. The decimal components and Majorana points are numerical renderings. Closed-form radical components are not claimed. A separate universal seventh-order exclusion proves maximum order six.',
        'verification':{'method':'Exact rational reduced-state equations, independently reconstructed Jacobian and strict Banach contraction bounds.',
            'exactExistenceMethod':'rational-contraction','exactNormalization':True,
            'exactOrderVerified':True,'globalMaximumOrderProven':True,'maximumOrderUpperBound':6,
            'existenceCertificate':'./research/n24-n26/n26-t6/certificate.json',
            'exclusionCertificate':'./research/n24-n26/n26-t7/certificate.json',
            'reductions':[{'t':q,'identityVerified':True,'purity':str(Fraction(1,q+1)),
                'evidence':'Certified rho_6=I_7/7, followed by partial trace.'} for q in range(1,7)]}}
    (HERE/'states.json').write_text(json.dumps({'states':[n24,n26]},indent=2,ensure_ascii=False)+'\n')
    print('Built two proven-existence records: radical N24 order7 and implicit-algebraic N26 order6.')

if __name__=='__main__':main()
