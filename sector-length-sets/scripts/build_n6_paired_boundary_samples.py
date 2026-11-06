#!/usr/bin/env python3
"""Explicit states on the sharp paired-face boundary; no optimization.

Existing scientific datasets remain read-only.  The real-Hodge construction is
from paired_face_sharp_upper_image_proof.txt; K(n) is the explicit seven-by-seven
matrix in N6_complete_paired_face_20261021.tex. Cofactors are evaluated at high
precision, and their nonzero chart is certified with interval arithmetic.
"""
from pathlib import Path
from fractions import Fraction as F
from math import sqrt, comb, dist
import json
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data'
AUDIT = ROOT/'audit'
mp.mp.dps = 90
mp.iv.dps = 45
LINEAR = [F(1,10),F(1,5),F(3,10),F(2,5),F(9,20),F(19,40)]
CURVED = [F(-1),F(-3,4),F(-1,2),F(-1,4),F(0),F(1,2),F(1),F(3,2),F(2),F(5,2)]
LINEAR_SOURCE = 'output/N6_threehour_border_search_20261020/side/paired_face_sharp_upper_image_proof.txt; output/pdf/assets/N6_complete_paired_face_20261021.tex'
CURVED_SOURCE = 'output/pdf/assets/N6_complete_paired_face_20261021.tex: sharp paired-face ellipse; direct odd-Dicke contraction'


def frac_tex(x):
    x = F(x)
    return str(x.numerator) if x.denominator == 1 else r'\frac{'+str(x.numerator)+'}{'+str(x.denominator)+'}'


def tuple6(x,y,z):
    return [x,y,z,45+10*x-2*y-3*z,3*z-9*x,18-2*x+y-z]


def linear_point(r):
    x=(1+4*r*r)/6
    return x,5+15*x,10-10*x


def ellipse_point(s):
    return (9-s*s)/24,(45+5*s)/4,(135-30*s-5*s*s)/24


def K_rows(n, ctx):
    a,b,c,d,e=n; z=ctx.mpf(0); q=ctx.sqrt
    return [
        [z,z,z,z,q(10)*e/5,-d,c],
        [z,z,z,-2*q(5)*e/5,q(15)*d/5,z,-b],
        [z,z,2*q(6)*e/5,-q(2)*d/5,-3*c/5,q(15)*b/5,q(10)*a/5],
        [z,-2*q(5)*e/5,-q(2)*d/5,4*c/5,-q(2)*b/5,-2*q(5)*a/5,z],
        [q(10)*e/5,q(15)*d/5,-3*c/5,-q(2)*b/5,2*q(6)*a/5,z,z],
        [-d,z,q(15)*b/5,-2*q(5)*a/5,z,z,z],
        [c,-b,q(10)*a/5,z,z,z,z]]


def kernel(r0, ctx):
    r=ctx.mpf(r0.numerator)/r0.denominator
    ell=(1+ctx.sqrt(1+4*r))/2
    h=r*(ell-1)
    U=ctx.mpf(1)/2-(ell-h)**2/3-h*h
    n0=(ell-h)/ctx.sqrt(6)
    n=[n0,(1/ctx.sqrt(2)+ctx.sqrt(U))/ctx.sqrt(2),h,
       (1/ctx.sqrt(2)-ctx.sqrt(U))/ctx.sqrt(2),n0]
    return n


def minor_rows(rows, row, col):
    return [[v for j,v in enumerate(values) if j!=col]
            for i,values in enumerate(rows) if i!=row]


def interval_det(rows):
    """Fixed-chart Gaussian interval determinant; every divisor excludes zero."""
    a=[list(row) for row in rows];det=mp.iv.mpf(1)
    for k in range(len(a)):
        possible=[i for i in range(k,len(a)) if a[i][k]>0 or a[i][k]<0]
        assert possible,'Uncertified interval pivot'
        pivot=max(possible,key=lambda i:abs(float(a[i][k].mid)))
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];det=-det
        p=a[k][k];det*=p
        for i in range(k+1,len(a)):
            ratio=a[i][k]/p
            for j in range(k+1,len(a)):a[i][j]-=ratio*a[k][j]
            a[i][k]=mp.iv.mpf(0)
    return det


def linear_state(r):
    n=kernel(r,mp.mp);rows=K_rows(n,mp.mp)
    B=[[rows[i][j]-(1 if i==j else 0) for j in range(7)] for i in range(7)]
    j=2
    cof=[(-1)**(j+k)*mp.det(mp.matrix(minor_rows(B,j,k))) for k in range(7)]
    scale=mp.sqrt(sum(c*c for c in cof));assert scale>mp.mpf('1e-5')
    c=[v/scale for v in cof]
    residual=max(abs(sum(rows[i][k]*c[k] for k in range(7))-c[i]) for i in range(7))
    assert residual<mp.mpf('1e-75')
    ni=kernel(r,mp.iv);ri=K_rows(ni,mp.iv)
    Bi=[[ri[i][k]-(1 if i==k else 0) for k in range(7)] for i in range(7)]
    chart=interval_det(minor_rows(Bi,j,j));assert chart>0
    return [float(v) for v in c],{'cofactorChart':j,'unitEigenvectorResidual':float(residual),
        'chartInterval':[float(chart.a),float(chart.b)],'intervalPrecisionDigits':mp.iv.dps}


def ellipse_state(s0):
    s=mp.mpf(s0.numerator)/s0.denominator
    A=mp.sqrt(3+mp.sqrt(9-s*s))
    return [0.,float(A/4),0.,float(mp.sqrt(10)/4),0.,float(s/(4*A)),0.]


def sectors(c):
    """All six sectors via Dicke-Schmidt marginal purities, independently."""
    out=[1.]
    for r in range(1,7):
        schmidt=[[c[k+j]*sqrt(comb(r,k)*comb(6-r,j)/comb(6,k+j))
            for j in range(7-r)] for k in range(r+1)]
        rho=[[sum(schmidt[k][j]*complex(schmidt[l][j]).conjugate() for j in range(7-r))
            for l in range(r+1)] for k in range(r+1)]
        purity=sum(abs(z)**2 for row in rho for z in row)
        out.append(comb(6,r)*(2**r*purity-sum(comb(r,k)/comb(6,k)*out[k] for k in range(r))))
    return [0. if abs(x)<1e-12 else x for x in out[1:]]


MATRIX_TEX = r'''\begin{pmatrix}
0&0&0&0&\sqrt{10}n_4/5&-n_3&n_2\\
0&0&0&-2\sqrt5n_4/5&\sqrt{15}n_3/5&0&-n_1\\
0&0&2\sqrt6n_4/5&-\sqrt2n_3/5&-3n_2/5&\sqrt{15}n_1/5&\sqrt{10}n_0/5\\
0&-2\sqrt5n_4/5&-\sqrt2n_3/5&4n_2/5&-\sqrt2n_1/5&-2\sqrt5n_0/5&0\\
\sqrt{10}n_4/5&\sqrt{15}n_3/5&-3n_2/5&-\sqrt2n_1/5&2\sqrt6n_0/5&0&0\\
-n_3&0&\sqrt{15}n_1/5&-2\sqrt5n_0/5&0&0&0\\
n_2&-n_1&\sqrt{10}n_0/5&0&0&0&0
\end{pmatrix}'''


def linear_expression(r):
    params=(r'r='+frac_tex(r)+r',\quad j=2\qquad\begin{aligned}'
        +r'\ell&=\frac{1+\sqrt{1+4r}}2,\quad h=r(\ell-1),\\'
        +r'U&=\frac12-\frac{(\ell-h)^2}{3}-h^2,\\'
        +r'n_0&=n_4=\frac{\ell-h}{\sqrt6},\quad n_2=h,\\'
        +r'n_1&=\frac{1/\sqrt2+\sqrt U}{\sqrt2},\quad '
        +r'n_3=\frac{1/\sqrt2-\sqrt U}{\sqrt2},\\'
        +r'K(n)&='+MATRIX_TEX+r',\\'
        +r'B&=K(n)-I_7,\\'
        +r'C_k&=(-1)^{j+k}\det B_{\widehat j,\widehat k},\quad 0\le k\le6,\\'
        +r'D_k&=|D_6^{(k)}\rangle.\end{aligned}')
    return {'latex':r'|L_r\rangle=\frac{\sum_{k=0}^{6}C_kD_k}{\sqrt{\sum_{k=0}^{6}C_k^2}}',
        'parameterLatex':params,
        'note':'Exact finite real-radical cofactor construction; no eigenvector fitting or interpolation.',
        'proofReference':LINEAR_SOURCE}


def ellipse_expression(s):
    return {'latex':r'|P_s\rangle=\frac{A_sD_1+\sqrt{10}D_3+(s/A_s)D_5}{4}',
        'parameterLatex':r's='+frac_tex(s)+r',\quad A_s=\sqrt{3+\sqrt{9-s^2}}',
        'note':'Exact normalized odd-C2 state on the sharp paired-face ellipse.',
        'proofReference':CURVED_SOURCE}


def atomic_text(path,text):
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(text);tmp.replace(path)


def main():
    originals=json.loads((DATA/'n6-manuscript-nodes.json').read_text())
    old=json.loads((DATA/'n6-boundary-samples.json').read_text())['samples']
    points=[p['point'] for p in originals]+[p['sectors'][:3] for p in old]
    samples=[];checks=[]
    for branch,params in [('linear',LINEAR),('ellipse',CURVED)]:
        for i,value in enumerate(params,1):
            if branch=='linear':
                c,chart=linear_state(value);target=linear_point(value);expression=linear_expression(value);symbol='L';parameter='r'
            else:
                c=ellipse_state(value);chart={};target=ellipse_point(value);expression=ellipse_expression(value);symbol='P';parameter='s'
            exact=tuple6(*target);actual=sectors(c)
            norm_error=abs(sum(x*x for x in c)-1)
            sector_error=max(abs(actual[k]-float(exact[k])) for k in range(6))
            nearest=min(dist([float(x) for x in target],p) for p in points)
            assert norm_error<2e-14 and sector_error<2e-10,(branch,value,norm_error,sector_error)
            assert nearest>1e-7,(branch,value,'duplicates a displayed node')
            sid=f'n6-paired-{branch}-{i}'
            samples.append({'id':sid,'N':6,'status':'attained','family':'Sharp paired-face boundary',
                'label':f'Paired-face boundary ({parameter} = {value})',
                'mathLabel':rf'$|{symbol}_{{{value}}}\rangle$',
                'parameters':{parameter:float(value)},'exactParameter':str(value),
                'coefficients':[[v,0.] for v in c],'sectors':actual,
                'proofOrigin':('Sharp paired-face boundary: real Hodge family' if branch=='linear' else 'Sharp paired-face boundary: odd C₂ ellipse'),
                'source':LINEAR_SOURCE if branch=='linear' else CURVED_SOURCE,
                'analyticalExpression':expression,
                'exactSectorValues':{'values':[frac_tex(v) for v in exact],'status':'exact',
                    'provenance':'Exact rational parameter in the proved sharp paired-face boundary image.'},
                'boundaryRole':'attained violet/gray boundary of the paired-purity face',
                'onOuterBoundary':True,'boundaryResidual':abs(15+5*actual[0]-actual[1]-actual[2])})
            checks.append({'id':sid,'parameter':str(value),'normError':norm_error,'sectorError':sector_error,
                'nearestPreviousNodeDistance':nearest,**chart})
    trace=[list(map(float,linear_point(F(j,240)))) for j in range(121)]
    trace += [list(map(float,ellipse_point(-1+F(j,180)))) for j in range(1,721)]
    assert len(trace)==841
    for sample in samples:
        assert min(dist(sample['sectors'][:3],p) for p in trace)<1e-10
    payload={'N':6,'samples':samples,'paths':[{'id':'paired-face-attained-boundary',
        'kind':'attained-state-locus','attained':True,'points':trace,
        'note':'Exact sharp paired-face boundary, from the chi_(1/2) contact through the odd ellipse to GHZ.'}],
        'notes':['All sixteen samples are actual normalized pure-state witnesses.',
            'The real-Hodge portion is a straight segment in sector coordinates; the odd-C2 portion is curved.',
            'No gray excluded outer crease is declared attainable. Existing nodes are unchanged.']}
    report={'status':'PASS','sampleCount':len(samples),'linearSamples':6,'ellipseSamples':10,'pathPoints':len(trace),
        'maximumNormError':max(c['normError'] for c in checks),
        'maximumSectorError':max(c['sectorError'] for c in checks),
        'allStatesOnTrace':True,'allPreviousNodesDistinct':True,'checks':checks,
        'proofNote':'audit/N6_paired_boundary_sample_proof.txt'}
    atomic_text(DATA/'n6-paired-boundary-samples.json',json.dumps(payload,indent=2,ensure_ascii=False)+'\n')
    atomic_text(DATA/'n6-paired-boundary-samples.js','globalThis.SECTOR_N6_PAIRED_BOUNDARY_SAMPLES='+json.dumps(payload,ensure_ascii=False,separators=(',',':'))+';\n')
    atomic_text(AUDIT/'n6-paired-boundary-samples-audit.json',json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))


if __name__=='__main__':main()
