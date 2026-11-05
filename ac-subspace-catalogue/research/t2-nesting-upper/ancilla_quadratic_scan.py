"""Exact inertia scan over irreducible ancillary spins and quadratic witnesses.
For each fixed N,s it checks every possible nonsingular sign pattern in the
one-parameter family X^2+aX-j(j+1)s(s+1)/3.
"""
from fractions import Fraction as F
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def best(n,twos):
    j=F(n,2);s=F(twos,2);c=j*(j+1)*s*(s+1)/3
    ls=[abs(j-s)+r for r in range(min(n,twos)+1)]
    eig=[(l*(l+1)-j*(j+1)-s*(s+1))/2 for l in ls]
    mult=[int(2*l+1)for l in ls];crit=sorted(set((c-x*x)/x for x in eig if x))
    trials=[crit[0]-1]+[(a+b)/2 for a,b in zip(crit,crit[1:])]+[crit[-1]+1]
    candidates=[]
    for a in trials:
        e=[x*x+a*x-c for x in eig];plus=sum(m for m,z in zip(mult,e)if z>0);minus=sum(m for m,z in zip(mult,e)if z<0)
        assert plus+minus==(n+1)*(twos+1)
        value=F(min(plus,minus),twos+1)
        candidates.append((value,a,plus,minus))
    return min(candidates)

def main():
    rows=[];improvements=[]
    for n in range(8,101):
        known=min(F(n-2,4),F(n+5,5));row=[]
        for twos in range(1,41):
            value,a,plus,minus=best(n,twos)
            row.append((value,twos,a,plus,minus))
            if value<known:improvements.append({'n':n,'ancillaTwiceSpin':twos,'bound':str(value),'a':str(a),'knownBound':str(known)})
        value,twos,a,plus,minus=min(row)
        rows.append({'n':n,'bestRealBound':str(value),'ancillaTwiceSpin':twos,'a':str(a),'inertia':[plus,minus]})
    data={'arithmetic':'Exact rational eigenvalues and all nonsingular sign intervals','maximumAncillaSpin':20,'improvements':improvements,'records':rows,'scope':'Finite scan within irreducible-ancilla rotationally invariant quadratic witnesses; not proof that no different general bound exists.'}
    (HERE/'ancilla-quadratic-scan.json').write_text(json.dumps(data,indent=2)+'\n')
    print('Improving rows',len(improvements));print(json.dumps(improvements[:30],indent=2))

if __name__=='__main__':main()
