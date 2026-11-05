#!/usr/bin/env python3
"""Exact Majorana multiplicities, including the projective point at infinity."""
import argparse,json,math,time
from pathlib import Path
from functools import reduce
import sympy as s
HERE=Path(__file__).resolve().parent
z=s.Symbol('z');y=s.Symbol('y');X,Y=s.symbols('X Y')

def analyze(state):
 started=time.time();n=state['n'];cs=[s.sympify(c) for c in state['exactCoefficients']];assert len(cs)==n+1
 coeff=[s.simplify((-1)**r*s.sqrt(s.binomial(n,r))*c) for r,c in enumerate(cs)]
 support=[r for r,c in enumerate(coeff) if c!=0];lo=min(support);hi=max(support);infinity=n-hi
 scale=coeff[hi];normalized=[s.simplify(c/scale) for c in coeff]
 step=reduce(math.gcd,[r-lo for r in support]) if len(support)>1 else 1
 h=s.Add(*[normalized[r]*y**((r-lo)//step) for r in support]);p=s.Poly(h,y,extension=True)
 content,factors=p.sqf_list();assert content==1
 reconstructed=s.Poly(1,y,domain=p.domain)
 factor_records=[];partition=([lo] if lo else [])+([infinity] if infinity else [])
 for f,m in factors:
  assert s.gcd(f,f.diff()).degree()==0
  for g,_ in factors:
   if f!=g:assert s.gcd(f,g).degree()==0
  reconstructed*=f**m
  partition.extend([int(m)]*(step*f.degree()))
  factor_records.append({'polynomial':s.sstr(f.as_expr()),'degree':f.degree(),'multiplicity':int(m),'nonzeroProjectiveRoots':step*f.degree()})
 assert reconstructed==p
 gcd=s.gcd(p,p.diff()) if p.degree()>0 else p
 finite=s.Add(*[normalized[r]*z**r for r in support]);hom=s.Add(*[normalized[r]*X**(n-r)*Y**r for r in support])
 partition.sort(reverse=True);assert sum(partition)==n
 record={'stateId':state['id'],'n':n,'achievedOrder':state['achievedOrder'],'label':state['label'],'partition':partition,'maxMultiplicity':max(partition),'distinctPoints':len(partition),'status':'exact','homogeneousPolynomial':s.sstr(hom),'homogeneousPolynomialLatex':s.latex(hom),'affinePolynomial':s.sstr(finite),'zeroMultiplicity':lo,'infinityMultiplicity':infinity,'nonzeroRootSubstitution':{'formula':f'z^{step}','step':step,'polynomial':s.sstr(p.as_expr()),'polynomialLatex':s.latex(p.as_expr())},'squarefreeFactors':factor_records,'gcdCertificate':{'polynomial':s.sstr(gcd.as_expr()),'degree':gcd.degree(),'domain':str(p.domain),'verifiedReconstruction':True,'verifiedSquarefreeFactors':True,'verifiedPairwiseCoprime':True},'references':state['references'],'runtimeSeconds':time.time()-started}
 print(state['id'],partition,'reduced',p.degree(),'gcd',gcd.degree(),'time',round(time.time()-started,3),flush=True)
 return record

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--input',type=Path,default=HERE/'states-input.json');parser.add_argument('--output',type=Path,default=HERE/'multiplicities.json');args=parser.parse_args();data=json.loads(args.input.read_text());records=[]
 for state in data['states']:
  records.append(analyze(state));args.output.write_text(json.dumps({'convention':{'homogeneous':'P(X,Y)=sum_r (-1)^r sqrt(binomial(n,r)) c_r X^(n-r)Y^r','affine':'z=Y/X; f(z)=P(1,z)','infinity':'[X:Y]=[0:1], multiplicity n-degree(f)','normalization':'Nonzero overall scalar omitted; leading affine coefficient set to1.','interpretation':'Multiplicity partition is independent of rotations, coordinate chart, complex conjugation and the antipodal convention.'},'method':'Exact algebraic squarefree factorization of f(z)=z^a h(z^g), including infinity. No numerical clustering.','states':records},indent=2))
if __name__=='__main__':main()
