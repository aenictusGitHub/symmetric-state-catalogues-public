"""Verify the full coherent/tangent Kramers family by exact polynomial substitution."""
import verify_certificate as v
from collections import defaultdict
from fractions import Fraction as F
P=defaultdict(F,{m:462*c for m,c in v.Q3.items()})
for m,c in v.product(v.Q,v.SK).items():P[m]-=c
R=defaultdict(F);mapping={0:(0,2),1:(1,2),10:(5,3),11:(4,3)}
for mon,c in P.items():
 ex=[0]*8
 for i in mon:
  if i%12 not in mapping:break
  if i%12==11:c=-c
  for j in mapping[i%12]:ex[(j+4*(i//12))%8]+=1
 else:R[tuple(ex)]+=c
assert not any(R.values())
print('PASS exact equality for all Kramers mixtures of the coherent/tangent plane.')
