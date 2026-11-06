"""Independent direct magnetic-index partial-time-reversal check.

This control compares a physical four-CG sum with the certificate's Racah
kernel. It is a convention audit, not the complete global certificate.
"""
import sys,json,random
from pathlib import Path
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
import flint
flint.ctx.prec=128
from angular import cg_ball
from validated_bootstrap import sixj_ball
A=flint.arb

def physical(n,S,T,J,L):
 v=A(0)
 for a in range(-n,n+1,2):
  ca=cg_ball(n,S,L,a,L-a)
  if ca.is_zero():continue
  for b in range(-n,n+1,2):
   cb=cg_ball(n,T,L,b,L-b)
   if cb.is_zero():continue
   phase=-1 if ((2*n-a-b)//2)%2 else 1
   v+=phase*ca*cb*cg_ball(n,S,J,-b,L-a)*cg_ball(n,T,J,-a,L-b)
 return A(flint.fmpq(L+1,J+1)).sqrt()*v

rng=random.Random(276);checks=[];n=27
for L in [0,2,4,6,8,10,12,14,16,18,20]:
 for k in range(3):
  S=rng.choice(list(range(max(1,n-L),n+L+1,2)))
  T=rng.choice(list(range(max(1,n-L),n+L+1,2)))
  J=rng.choice(list(range(max(abs(n-S),abs(n-T)),min(n+S,n+T)+1,2)))
  v=physical(n,S,T,J,L)
  target=A((J+1)*(L+1)).sqrt()*sixj_ball(n,S,J,n,T,L)
  diff=v-target
  assert diff.abs_upper()<A('1e-28'),(S,T,J,L,str(diff))
  checks.append(dict(S=S,T=T,J=J,L=L,errorUpperRational=str(diff.abs_upper().fmpq())))
r=dict(passed=True,n=n,copies=6,split=[1,5],precisionBits=128,checks=checks,scope='Direct magnetic-index four-CG comparison; phase and normalization audit only, separate from the global validated certificate.')
(P/'recoupling-control.json').write_text(json.dumps(r,indent=2)+'\n')
print('PASS',len(checks),'N27 checks across every used output spin')
