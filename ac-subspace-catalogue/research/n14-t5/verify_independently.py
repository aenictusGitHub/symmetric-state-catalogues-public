#!/usr/bin/env python3
"""Independent integer audit of a selected three-copy AC exclusion certificate.

Uses symmetrized tensor words and explicit single-word angular momentum;
imports no discovery code. Certificate schema follows the N11 moment proof.
"""
from pathlib import Path
from itertools import combinations_with_replacement
from collections import Counter
from math import comb,gcd,factorial,isqrt,prod
import json,argparse,hashlib
from fractions import Fraction
from functools import lru_cache
HERE=Path(__file__).resolve().parent

def modular_rank(rows,p):
    a=[[int(x)%p for x in row] for row in rows];rank=0;det=1
    for col in range(len(a[0])):
        pivot=next((i for i in range(rank,len(a)) if a[i][col]),None)
        if pivot is None:continue
        if pivot!=rank:a[pivot],a[rank]=a[rank],a[pivot];det=-det
        value=a[rank][col];det=det*value%p;inv=pow(value,-1,p)
        a[rank]=[(z*inv)%p for z in a[rank]]
        for i in range(rank+1,len(a)):
            factor=a[i][col]
            if factor:a[i]=[(x-factor*y)%p for x,y in zip(a[i],a[rank])]
        rank+=1
        if rank==len(a):break
    return rank,det%p

@lru_cache(None)
def polynomials(n,t):
    values=[]
    for ell in range(t+1):
        v=[Fraction((2*i-n)**ell)for i in range(n+1)]
        for u in values:
            ratio=sum(a*b for a,b in zip(v,u))/sum(a*a for a in u)
            v=[a-ratio*b for a,b in zip(v,u)]
        den=1
        for x in v:den=den*x.denominator//gcd(den,x.denominator)
        vv=[int(x*den)for x in v];common=gcd(*vv);values.append(tuple(x//common for x in vv))
    return tuple(values)

def f(l,i,n):return polynomials(n,5)[l][i]

def verify(path=None):
    path=Path(path or HERE/'certificate.json');raw=path.read_bytes();cert=json.loads(raw)
    n,r,t=cert['n'],cert['copies'],cert['t'];assert(n,r,t)==(14,3,5)
    assert {'rank','kernel','positiveDiagonalCertificate'}<=set(cert)
    words=list(combinations_with_replacement(range(n+1),r))
    byweight={}
    for word in words:byweight.setdefault(2*sum(word)-r*n,[]).append(word)
    offsets={w:{word:i for i,word in enumerate(group)}for w,group in byweight.items()}
    lower=list(combinations_with_replacement(range(n+1),r-1));lowweight={}
    for word in lower:lowweight.setdefault(2*sum(word)-(r-1)*n,[]).append(word)
    prime=int(cert.get('determinantModulus',1000003))
    assert prime>2 and all(prime%p for p in range(2,isqrt(prime)+1))
    spins=[int(record['J2'])for record in cert['blocks']]
    assert spins==sorted(set(spins)) and all(J>=0 for J in spins)
    blocks=[];total=0
    for record in cert['blocks']:
        J=int(record['J2']);mult=int(record['multiplicity'])
        highest=[[int(x)for x in row]for row in record['highestWeight']]
        assert len(highest)==len(byweight[J]) and all(len(row)==mult for row in highest)
        assert modular_rank(highest,prime)[0]==mult
        def act(matrix,M,direction):
            out=[[0]*mult for _ in byweight.get(M+2*direction,[])]
            for word,row in zip(byweight[M],matrix):
                for i,count in Counter(word).items():
                    nxt=i+direction
                    if not 0<=nxt<=n:continue
                    replacement=list(word);replacement.remove(i);replacement.append(nxt);replacement=tuple(sorted(replacement))
                    index=offsets[M+2*direction][replacement]
                    factor=count*((i+1)if direction==1 else(n-i+1))
                    for a,value in enumerate(row):out[index][a]+=factor*value
            return out
        assert all(not any(row)for row in act(highest,J,1))
        descendants={J:highest}
        for step in range(1,J+1):
            M=J-2*step;v=act(descendants[M+2],M+2,-1)
            assert all(x%step==0 for row in v for x in row)
            v=[[x//step for x in row]for row in v];descendants[M]=v
            assert act(v,M,1)==[[x*(J-step+1)for x in row]for row in descendants[M+2]]
        den=int(record['denominatorScale']);assert den>0 and all(den%comb(J,k)==0 for k in range(J+1))
        blocks.append((J,mult,den,descendants));total+=(J+1)*mult
    assert total==len(words)==680
    variables=[(ib,a,b)for ib,(_,m,_,_)in enumerate(blocks)for a in range(m)for b in range(a,m)]
    assert variables==[tuple(v)for v in cert['variables']]
    for l in range(t+1):
        for h in range(l):assert sum(f(l,i,n)*f(h,i,n)for i in range(n+1))==0
    rebuilt=[]
    for ell,M,ia,ib in cert['selectedConstraints']:
        assert 1<=ell<=t and 0<=ia<=ib<len(lowweight[M])
        a,b=lowweight[M][ia],lowweight[M][ib];ca,cb=Counter(a),Counter(b)
        # Exactly the first tensor factor's contraction, using orbit multiplicities.
        terms=[]
        for i in range(n+1):
            weight=M+2*i-n;aa=tuple(sorted(a+(i,)));bb=tuple(sorted(b+(i,)))
            factor=comb(n,i)*(ca[i]+1)*(cb[i]+1)*f(ell,i,n)
            terms.append((weight,offsets[weight][aa],offsets[weight][bb],factor))
        row=[]
        for blockindex,alpha,beta in variables:
            J,_,den,U=blocks[blockindex];value=0
            for weight,aidx,bidx,factor in terms:
                if weight not in U:continue
                aa,bb=U[weight][aidx],U[weight][bidx]
                product=aa[alpha]*bb[beta]
                if alpha!=beta:product+=aa[beta]*bb[alpha]
                value+=den//comb(J,(J-weight)//2)*factor*product
            row.append(value)
        common=gcd(*row);assert common>0;rebuilt.append([x//common for x in row])
    assert rebuilt==[[int(x)for x in row]for row in cert['selectedRows']]
    rank,residue=modular_rank(rebuilt,prime)
    assert rank==len(rebuilt)==cert.get('rank',len(variables))
    if 'kernel' in cert:
        kernel=[[int(x)for x in row]for row in cert['kernel']]
        assert len(kernel)==len(variables) and all(len(row)==len(variables)-rank for row in kernel)
        assert modular_rank(kernel,prime)[0]==len(variables)-rank
        assert all(sum(a*b for a,b in zip(row,column))==0 for row in rebuilt for column in zip(*kernel))
    else:
        assert rank==len(variables) and residue!=0
        assert residue==cert['determinantResidue']
    positive=cert['positiveDiagonalCertificate'];indices=positive['indices']
    primitive=[]
    for index in indices:
        _,a,b=variables[index];assert a==b
        row=kernel[index];common=gcd(*row);assert common>0
        primitive.append([x//common for x in row])
    assert primitive==positive['primitiveRows']
    weights=[int(x)for x in positive['positiveWeights']]
    assert len(weights)==len(primitive)==3 and all(x>0 for x in weights)
    assert all(sum(w*row[col]for w,row in zip(weights,primitive))==0 for col in range(2))
    assert modular_rank(primitive,prime)[0]==2
    # Independent physical probabilities, reconstructed using every descendant.
    prob=cert['probabilityCertificate'];probrows=[]
    assert prob['spins']==[4,14,16] and len(prob['blockIndices'])==3
    metrics={M:[prod(comb(n,i)**count*factorial(count)for i,count in Counter(word).items())
        for word in group]for M,group in byweight.items()}
    for spin,blockindex in zip(prob['spins'],prob['blockIndices']):
        J,mult,den,U=blocks[blockindex];assert J==2*spin
        block=[[0]*mult for _ in range(mult)]
        for M,matrix in U.items():
            scale=den//comb(J,(J-M)//2)
            for metric,row in zip(metrics[M],matrix):
                for a in range(mult):
                    for b in range(mult):block[a][b]+=scale*metric*row[a]*row[b]
        H=U[J]
        highest=[[den*(J+1)*sum(metric*row[a]*row[b]for metric,row in zip(metrics[J],H))
            for b in range(mult)]for a in range(mult)]
        assert block==highest
        form=[Fraction(0),Fraction(0)]
        for index,(ib,a,b)in enumerate(variables):
            if ib!=blockindex:continue
            coefficient=Fraction(block[a][b]*(1 if a==b else 2),factorial(r))
            for col in range(2):form[col]+=coefficient*kernel[index][col]
        probrows.append(form)
    assert probrows==[[Fraction(x)for x in row]for row in prob['kernelRows']]
    weights_p=[int(x)for x in prob['positiveWeights']];assert all(x>0 for x in weights_p)
    assert all(sum(w*row[col]for w,row in zip(weights_p,probrows))==0 for col in range(2))
    primitive_p=[]
    for row in probrows:
        den=1
        for x in row:den=den*x.denominator//gcd(den,x.denominator)
        ir=[int(x*den)for x in row];common=gcd(*ir);assert common>0
        primitive_p.append([x//common for x in ir])
    assert modular_rank(primitive_p,prime)[0]==2
    # A direct trace-versus-probabilities contradiction, without rank reasoning.
    totaltrace=[Fraction(0),Fraction(0)]
    for ib,(J,mult,den,U)in enumerate(blocks):
        block=[[0]*mult for _ in range(mult)]
        for M,matrix in U.items():
            scale=den//comb(J,(J-M)//2)
            for metric,row in zip(metrics[M],matrix):
                for a in range(mult):
                    for b in range(mult):block[a][b]+=scale*metric*row[a]*row[b]
        H=U[J]
        assert block==[[den*(J+1)*sum(metric*row[a]*row[b]for metric,row in zip(metrics[J],H))
            for b in range(mult)]for a in range(mult)]
        for index,(jb,a,b)in enumerate(variables):
            if jb!=ib:continue
            coefficient=Fraction(block[a][b]*(1 if a==b else 2),factorial(r))
            for col in range(2):totaltrace[col]+=coefficient*kernel[index][col]
    traceprob=cert['traceProbabilityCertificate']
    assert totaltrace==[Fraction(x)for x in traceprob['traceKernelRow']]
    tcoef=int(traceprob['traceCoefficient']);pc=traceprob['positiveProbabilityCoefficients']
    assert tcoef>0 and all(x>0 for x in pc) and traceprob['spins']==[4,16]
    selected=[probrows[prob['spins'].index(j)]for j in traceprob['spins']]
    assert all(tcoef*totaltrace[col]+sum(c*row[col]for c,row in zip(pc,selected))==0 for col in range(2))
    result={'passed':True,'method':'Independent integer symmetric-word reconstruction; no discovery imports.',
        'n':n,'t':t,'copies':r,'momentDimension':len(words),'realInvariantVariables':len(variables),
        'constraintRank':rank,'modularRankPrime':prime,'kernelDimension':len(variables)-rank,
        'scope':'Universal exclusion for arbitrary complex pure symmetric states.',
        'certificateSHA256':hashlib.sha256(raw).hexdigest(),
        'positiveDiagonalCertificate':{'indices':indices,'primitiveRows':primitive,'strictlyPositiveWeights':weights,
            'positiveCombinationVanishesExactly':True,'formsHaveRankTwo':True,
            'conclusion':'Positivity forces each of the three diagonals to vanish; rank two forces the whole moment to vanish, contradicting trace one.'},
        'probabilityCertificate':{'spins':prob['spins'],'strictlyPositiveWeights':weights_p,
            'probabilitiesReconstructedFromAllWords':True,'highestWeightShortcutIndependentlyChecked':True,
            'positiveCombinationVanishesExactly':True,'formsHaveRankTwo':True,'identity':prob['identity']},
        'traceProbabilityCertificate':{'identity':traceprob['identity'],
            'physicalTraceReconstructedFromAllWords':True,'highestWeightShortcutIndependentlyChecked':True,
            'strictlyPositiveCoefficients':[tcoef]+pc,'identityVerifiedExactly':True,
            'contradiction':'Trace one and nonnegative probabilities make the left side strictly positive.'},
        'statement':'No pure symmetric fourteen-qubit state is fifth-order anticoherent: K(14,5)=0.'}
    (HERE/'independent-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('certificate',nargs='?');args=p.parse_args();print(json.dumps(verify(args.certificate),indent=2))
