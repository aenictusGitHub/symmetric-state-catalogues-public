"""Verify the N=11, t=4 stationary state and its strict local minimum modulo rotations."""
import json
import sys
from math import comb
from fractions import Fraction
from pathlib import Path

from sage.all import ComplexIntervalField, QQ, RDF, RealIntervalField, identity_matrix, matrix, vector

R = RealIntervalField(160)
C = ComplexIntervalField(160)
I = C.gen()
root = Path(__file__).resolve().parent
z_center = json.loads((root / 'witness_t4.json').read_text())['z']
even = [(2, 0)] + [(k, part) for k in range(4, 11, 2) for part in (0, 1)]
odd = [(k, part) for k in range(3, 12, 2) for part in (0, 1)]
center = [QQ(Fraction(z_center[k][part])) for k, part in even]
radius = QQ(1) / 10**36
box = [R(v-radius, v+radius) for v in center]
weights = matrix(R, 5, 8, lambda i, j: R(comb(4, i)*comb(7, j)).sqrt())

def adjoint(A):
    return A.conjugate().transpose()

def bipartition(z):
    return matrix(C, 5, 8, lambda i, j: weights[i, j]*z[i+j])

def direction(k, part):
    z = [C(0)]*12
    z[k] = I if part else C(1)
    return bipartition(z)

def derivatives(coords, chart):
    z = [C(0)]*12
    z[0] = C(1)
    for (k, part), value in zip(even, coords):
        z[k] += I*value if part else value
    A = bipartition(z)
    rho = A*adjoint(A)
    n, F = rho.trace().real(), (rho*rho).trace().real()
    H = [direction(*axis) for axis in chart]
    B = [h*adjoint(A)+A*adjoint(h) for h in H]
    n1 = [b.trace().real() for b in B]
    F1 = [2*(rho*b).trace().real() for b in B]
    n2, F2 = matrix(R, len(H)), matrix(R, len(H))
    for i, h in enumerate(H):
        for j in range(i+1):
            p = h*adjoint(H[j])
            cross = p+adjoint(p)
            n2[i,j] = n2[j,i] = cross.trace().real()
            F2[i,j] = F2[j,i] = 2*(B[i]*B[j]+rho*cross).trace().real()
    g = vector(R, [n*f-2*F*d for f,d in zip(F1,n1)])
    jac = matrix(R, len(H), len(H), lambda i,j:
                 n*F2[i,j]+n1[j]*F1[i]-2*F1[j]*n1[i]-2*F*n2[i,j])
    beta = F/n**2
    hess = matrix(R, len(H), len(H), lambda i,j:
                  F2[i,j]-2*beta*(n1[i]*n1[j]+n*n2[i,j]))
    return beta, g, jac, hess

_, g0, J0, _ = derivatives(center, even)
beta, _, J, _ = derivatives(box, even)
m = len(even)
approx = matrix(RDF, m, m, lambda i,j: float(J0[i,j].center())).inverse()
preconditioner = approx.change_ring(R)
def norm_inf(A):
    return max(sum(abs(v) for v in row).upper() for row in A.rows())

contraction = norm_inf(identity_matrix(R, m)-preconditioner*J)
defect = max(abs(v).upper() for v in preconditioner*g0)
assert (R(defect)+R(contraction)*R(radius)).upper() < R(radius).lower()
# The odd first derivatives vanish by parity. The transverse rotation
# generators act on z1 with real-linear factors proportional to 1 ± 10*z2.
assert box[0].lower() > 0 and not (1-10*box[0]).contains_zero()
_, _, _, hess = derivatives(box, even+odd)
minimum_pivot = R(1).lower()
for k in range(hess.nrows()):
    pivot = hess[k,k]
    assert pivot.lower() > 0
    minimum_pivot = min(minimum_pivot, pivot.lower())
    for i in range(k+1, hess.nrows()):
        for j in range(k+1, hess.nrows()):
            hess[i,j] -= hess[i,k]*hess[k,j]/pivot

print('Contraction bound:', contraction)
print('Minimum Hessian pivot lower bound:', minimum_pivot)
print('Stationary purity enclosure:', beta)
print('Verified a strict local minimum modulo collective rotations and phase.')

if '--json' in sys.argv:
    # Preserve the tighter rational endpoints used by the interpolation run.
    # QQ(MPFR) may round to a simple rational. Verify the exported box against
    # a separately enclosed smaller fixed-point box before using it.
    inner_radius = radius/2
    assert (R(defect)+R(contraction)*R(inner_radius)).upper() < R(inner_radius).lower()
    inner_box = [R(v-inner_radius,v+inner_radius) for v in center]
    inner_beta, *_ = derivatives(inner_box,even)
    def certified_bounds(outer, inner):
        lo, hi = QQ(outer.lower()), QQ(outer.upper())
        assert lo <= inner.lower().exact_rational()
        assert inner.upper().exact_rational() <= hi
        return [str(lo),str(hi)]
    print(json.dumps(dict(even=even,
          box=[certified_bounds(a,b) for a,b in zip(box,inner_box)],
          beta=certified_bounds(beta,inner_beta))))
