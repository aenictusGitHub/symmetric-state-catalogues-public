#!/usr/bin/env python3
"""Verify every exact AC witness in states.json; requires SymPy only."""
import argparse,json
from pathlib import Path
import sympy as s

def verify(path):
 data=json.loads(path.read_text())
 for entry in data['states']:
  n=entry['n'];order=entry['achievedOrder']
  c=[s.sympify(x,locals={'I':s.I,'sqrt':s.sqrt}) for x in entry['exactCoefficients']]
  assert len(c)==n+1,entry['id']
  assert s.simplify(sum(x*s.conjugate(x) for x in c))==1,entry['id']
  for t in range(1,min(order+1,n//2)+1):
   M=s.Matrix(t+1,n-t+1,lambda a,b:c[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(n-t,b)/s.binomial(n,a+b)))
   delta=(M*M.conjugate().T-s.eye(t+1)/s.Integer(t+1)).applyfunc(s.simplify)
   is_zero=all(x==0 for x in delta)
   assert is_zero==(t<=order),(entry['id'],t,delta)
  print(f"PASS {entry['id']}: normalized; full AC order {order}; all lower orders verified")
 print(f"Verified {len(data['states'])} exact witnesses. No claim about maximal possible order is inferred.")

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('path',nargs='?',type=Path,default=Path(__file__).resolve().parent/'states.json')
 verify(parser.parse_args().path)
