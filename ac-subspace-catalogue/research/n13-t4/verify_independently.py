#!/usr/bin/env python3
"""Independent integer audit of a selected four-copy AC exclusion certificate.

Uses symmetrized tensor words and explicit single-word angular momentum;
imports no discovery code. Certificate schema follows the N11 moment proof.
"""
from pathlib import Path
from itertools import combinations_with_replacement
from collections import Counter
from math import comb,gcd,factorial,isqrt,prod
import json,argparse,hashlib
from fractions import Fraction
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

def f(l,i,n):
    z=2*i-n;v=n*(n+2)
    return [1,z,3*z*z-v,z*(5*z*z-3*v+4),35*z**4-(30*v-100)*z*z+3*v*(v-8)][l]

def verify(path=None):
    path=Path(path or HERE/'certificate.json');raw=path.read_bytes();cert=json.loads(raw)
    n,r,t=cert['n'],cert['copies'],cert['t'];assert(n,r,t)==(13,4,4)
    assert {'rank','kernel','traceRow','traceBlocks'}<=set(cert)
    words=list(combinations_with_replacement(range(n+1),r))
    byweight={}
    for word in words:byweight.setdefault(2*sum(word)-r*n,[]).append(word)
    offsets={w:{word:i for i,word in enumerate(group)}for w,group in byweight.items()}
    lower=list(combinations_with_replacement(range(n+1),r-1));lowweight={}
    for word in lower:lowweight.setdefault(2*sum(word)-(r-1)*n,[]).append(word)
    prime=int(cert.get('determinantModulus',1000003))
    assert prime>2 and all(prime%p for p in range(2,isqrt(prime)+1))
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
    assert total==len(words)==2380
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
    trace_check=None
    if 'traceRow' in cert:
        # Reconstruct the trace from every physical symmetrized word.  A word
        # orbit has metric prod(binomial(n,i)^occ_i * occ_i!)/r!.
        metrics={M:[prod(comb(n,i)**count*factorial(count)for i,count in Counter(word).items())
            for word in group] for M,group in byweight.items()}
        traceblocks=[];tracecoeffs=[]
        for J,mult,den,U in blocks:
            block=[[0]*mult for _ in range(mult)]
            for M,statesmatrix in U.items():
                scale=den//comb(J,(J-M)//2)
                for metric,row in zip(metrics[M],statesmatrix):
                    for a in range(mult):
                        for b in range(mult):block[a][b]+=scale*metric*row[a]*row[b]
            # Independently check the highest-weight shortcut in the primary proof.
            H=U[J]
            highest=[[den*(J+1)*sum(metric*row[a]*row[b]for metric,row in zip(metrics[J],H))
                for b in range(mult)]for a in range(mult)]
            assert block==highest
            block=[[Fraction(value,factorial(r))for value in row]for row in block]
            traceblocks.append(block)
            tracecoeffs.extend(block[a][b]*(1 if a==b else 2)for a in range(mult)for b in range(a,mult))
        assert traceblocks==[[[Fraction(x)for x in row]for row in block]for block in cert['traceBlocks']]
        assert tracecoeffs==[Fraction(x)for x in cert['traceRow']]
        assert all(sum(c*x for c,x in zip(tracecoeffs,column))==0 for column in zip(*kernel))
        trace_check={'physicalTraceReconstructedFromAllWords':True,'highestWeightShortcutIndependentlyChecked':True,
            'traceAnnihilatesWholeKernelExactly':True,'contradiction':'The invariant fourth moment of a normalized pure state has trace one; the exact constraints force trace zero.'}
    result={'passed':True,'method':'Independent integer symmetric-word reconstruction; no discovery imports.',
        'n':n,'t':t,'copies':r,'momentDimension':len(words),'realInvariantVariables':len(variables),
        'constraintRank':rank,'modularRankPrime':prime,'kernelDimension':len(variables)-rank,
        'scope':'Necessary invariant moment constraints only; nonzero kernel needs a separate positivity exclusion.',
        'certificateSHA256':hashlib.sha256(raw).hexdigest()}
    if trace_check:
        result['traceCertificate']=trace_check
        result['statement']='No pure symmetric thirteen-qubit state is fourth-order anticoherent: K(13,4)=0.'
        result['scope']='Universal exclusion for arbitrary complex pure symmetric states, using an exact trace identity.'
    (HERE/'independent-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('certificate',nargs='?');args=p.parse_args();print(json.dumps(verify(args.certificate),indent=2))
