#!/usr/bin/env python3
"""Verify full-subspace anticoherence, including all off-diagonal cross terms.
Uses only normalized Dicke basis splitting, independently of tensor conventions.
"""
from pathlib import Path
import json, math
import sympy as s
import numpy as np
ROOT=Path(__file__).resolve().parent

def split_matrix(c,n,t):
    return s.Matrix(t+1,n-t+1,lambda a,b:s.simplify(c[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(n-t,b)/s.binomial(n,a+b))))

def symbolic_check(rows,n,t):
    vectors=[s.Matrix([s.sympify(v) for v in row]) for row in rows]
    gram=s.simplify(s.Matrix.hstack(*vectors).H*s.Matrix.hstack(*vectors));orthonormal=gram==s.eye(len(rows))
    matrices=[split_matrix(c,n,t) for c in vectors]
    bad=[];checked=0
    for a,A in enumerate(matrices):
        for b,B in enumerate(matrices):
            deviation=A*B.H-gram[b,a]*s.eye(t+1)/(t+1)
            for u in range(t+1):
                for v in range(t+1):
                    checked+=1;value=s.simplify(s.expand_complex(deviation[u,v]))
                    if value!=0:
                        bad.append({'basisA':a,'basisB':b,'reducedRow':u,'reducedColumn':v,'value':str(value)})
    return {'passed':not bad,'orthonormal':orthonormal,'fullRank':gram.det()!=0,'matrixElementsChecked':checked,'violations':bad[:4]}

def numeric_check(rows,n,t):
    B=np.array([[complex(*v) for v in row] for row in rows],complex).T
    gram=B.conj().T@B
    matrices=[np.array([[B[a+b,c]*math.sqrt(math.comb(t,a)*math.comb(n-t,b)/math.comb(n,a+b)) for b in range(n-t+1)] for a in range(t+1)]) for c in range(B.shape[1])]
    deviations=[A@D.conj().T-gram[d,a]*np.eye(t+1)/(t+1) for a,A in enumerate(matrices) for d,D in enumerate(matrices)]
    return {'gramError':float(np.linalg.norm(gram-np.eye(B.shape[1]))),'maxCrossReductionFrobeniusError':max(float(np.linalg.norm(x)) for x in deviations),'combinedCrossReductionFrobeniusError':math.sqrt(sum(float(np.linalg.norm(x))**2 for x in deviations)),'allOrderedBasisPairsChecked':len(deviations),'rank':int(np.linalg.matrix_rank(B))}

def main():
 data=json.loads((ROOT/'findings.json').read_text());reports=[]
 for r in data['exactWitnesses']:
  check=symbolic_check(r['exactBasis'],r['n'],r['order']);assert check['passed'] and check['orthonormal'] and check['fullRank'],(r['id'],check)
  achieved=r['order'];higher=symbolic_check(r['exactBasis'],r['n'],achieved+1)
  while higher['passed'] and achieved<r['n']//2:
   achieved+=1;higher=symbolic_check(r['exactBasis'],r['n'],achieved+1)
  r['verifiedAchievedOrder']=achieved
  row={'id':r['id'],'n':r['n'],'dimension':r['dimension'],'order':r['order'],'verifiedAchievedOrder':achieved,**check,'nextOrderFails':not higher['passed'],'nextOrderCounterexample':higher['violations'][0] if higher['violations'] else None}
  reports.append(row);print(r['id'],'verified',flush=True)
 numerical=[{'id':r['id'],**numeric_check(r['basis'],r['n'],r['order'])} for r in data['numericalWitnesses']]
 # A concrete counterexample to adding the center vector at odd integer j.
 wrong=[['1/sqrt(2)' if r in [0,6] else '0' for r in range(7)],['1/sqrt(2)' if r in [2,4] else '0' for r in range(7)],['1' if r==3 else '0' for r in range(7)]]
 caveat=symbolic_check(wrong,6,1);assert not caveat['passed']
 out={'passed':True,'method':'Exact Dicke splitting matrices A_a and every ordered cross product A_a A_b† = <ψ_b|ψ_a>I/(t+1), including all basis pairs. Normalization, orthogonality and full rank checked. Next order tested separately.','exactWitnesses':reports,'numericalWitnesses':numerical,'eq48OddIntegerCounterexample':caveat,'limitations':['Numerical residuals do not certify exact existence.','An attained dimension is a lower bound until it meets an independent upper bound.']}
 (ROOT/'findings.json').write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
 (ROOT/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'exactWitnesses':len(reports),'exactMatrixElementsChecked':sum(r['matrixElementsChecked'] for r in reports),'numericalWitnesses':numerical,'allClaimedOrdersVerified':True},indent=2))
if __name__=='__main__':main()
