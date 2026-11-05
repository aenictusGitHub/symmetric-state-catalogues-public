"""Build certified N23 and radical N25 catalogue records; mpmath for rendering."""
from pathlib import Path
from fractions import Fraction
import copy,json
import mpmath as mp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def main():
    source=json.loads((ROOT/'data/catalogue.json').read_text())
    original=next(s for s in source['states'] if s['id']=='liege-acs.23.5')
    cert=json.loads((HERE/'n23-t5/certificate.json').read_text())
    center=list(map(Fraction,cert['center']));scales=cert['scaleFactors'];mp.mp.dps=150
    rendered=[]
    for r in range(24):
        factor=mp.sqrt(mp.binomial(23,r))/scales[r]
        rendered.append([mp.nstr(factor*mp.mpf(center[k].numerator)/center[k].denominator,100)
            for k in [r,24+r]])
    n23={'id':'certified-algebraic-n23-t5','n':23,'achievedOrder':5,'witnessesOrders':list(range(1,6)),
        'label':'Certified algebraic state · order 5','status':'proven-existence',
        'coefficientEvidence':'implicit-algebraic',
        'stateLatex':r'|\psi\rangle=\sum_{r=0}^{23}\frac{\sqrt{\binom{23}{r}}}{d_r}\bigl(x_r^*+ix_{24+r}^*\bigr)|D_{23}^{(r)}\rangle',
        'support':list(range(24)),'coefficients':[[float(x),float(y)] for x,y in rendered],
        'highPrecisionCoefficients':rendered,
        'majoranaPoints':copy.deepcopy(original['majoranaPoints']),'partition':copy.deepcopy(original['partition']),
        'majoranaCoordinatePolicy':'Numeric rendering copied from the nearby refined Liège state; the exact state is defined by the certified polynomial root.',
        'normalization':'Exact rho_5 = I_6/6 implies sum_r |c_r|^2 = 1.',
        'coefficientConvention':'Normalized Dicke basis; r counts excitations. Floating components render the rational certificate center.',
        'implicitAlgebraicDefinition':{
            'certificateUrl':'./research/n23-n25/n23-t5/certificate.json','proofUrl':'./research/n23-n25/n23-t5/proof.tex',
            'equations':'36 real rational quadratics equivalent to rho_5=I_6/6 after the stated binomial scaling.',
            'variableConvention':cert['variableConvention'],'activeVariables':cert['activeVariables'],
            'fixedVariables':cert['fixedVariables'],'radius':cert['radius'],
            'rootDefinition':'Unique real zero in the active-coordinate infinity-norm box around the saved rational center, with twelve coordinates fixed to their saved rational values.',
            'centerCoefficientErrorBound':'2e-65 per real or imaginary Dicke component, relative to the exact rational-center conversion; ordinary floating values have additional rounding.',
            'closedFormRadicalsClaimed':False},
        'derivedFromStateId':'liege-acs.23.5','sourceLibraryRecordId':'acs.23.5',
        'references':[{'title':'Twenty-three qubits: certified maximum order five','url':'./research/n23-n25/index.html#n23'},
            {'title':'Exact algebraic existence proof','url':'./research/n23-n25/n23-t5/proof.tex'},
            {'title':'Rational root certificate','url':'./research/n23-n25/n23-t5/certificate.json'},
            {'title':'Independent rational verification','url':'./research/n23-n25/n23-t5/independent-verification.json'},
            {'title':'Universal sixth-order exclusion','url':'./research/n23-n25/n23-t6/proof.tex'}],
        'provenance':'Certified exact algebraic zero near refined Liège acs.23.5; exact rational Banach contraction certificate.',
        'notes':'Exact fifth-order existence is proved by an isolating rational contraction certificate. Decimal components and Majorana points are numerical renderings. Closed-form radical components are not claimed. A separate universal sixth-order exclusion proves maximum order five.',
        'verification':{'method':'Exact rational reduced-state equations, independently reconstructed Jacobian and strict Banach contraction bounds.',
            'exactExistenceMethod':'rational-contraction','exactNormalization':True,'exactOrderVerified':True,
            'globalMaximumOrderProven':True,'maximumOrderUpperBound':5,
            'existenceCertificate':'./research/n23-n25/n23-t5/certificate.json',
            'exclusionCertificate':'./research/n23-n25/n23-t6/certificate.json',
            'reductions':[{'t':q,'identityVerified':True,'purity':str(Fraction(1,q+1)),
                'evidence':'Certified rho_5=I_6/6, followed by partial trace.'} for q in range(1,6)]}}
    n25=copy.deepcopy(json.loads((HERE/'n25-t5/states.json').read_text())['states'][0])
    n25['stateLatex']=r'\begin{gathered}'+n25['stateLatex']+r'\\[6pt]'+n25['parametersLatex']+r'\end{gathered}'
    n25['majoranaCoordinatePolicy']='numeric-roots-retain-explicit-poles'
    n25['references']=[{'title':'Twenty-five qubits: exact order five and rigorous bounds','url':'./research/n23-n25/index.html#n25'},
        {'title':'Analytical radical construction','url':'./research/n23-n25/n25-t5/proof.tex'},
        {'title':'Independent exact reduced-state verification','url':'./research/n23-n25/n25-t5/independent-verification.json'},
        {'title':'Universal seventh-order exclusion','url':'./research/n23-n25/n25-t7/proof.tex'}]
    n25['verification']['maximumOrderUpperBound']=6
    n25['notes']='Exactly fifth-order AC. A universal seventh-order exclusion gives maximum order at most six. Sixth-order existence at N=25 remains unresolved; no maximum-five claim is made.'
    (HERE/'states.json').write_text(json.dumps({'states':[n23,n25]},indent=2,ensure_ascii=False)+'\n')
    print('Built certified maximum-five N23 and exact fifth-order N25 records; N25 remains 5–6.')

if __name__=='__main__':main()
