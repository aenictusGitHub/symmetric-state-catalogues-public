"""Exact rational verification of the N=12 half-system purity certificate."""

from fractions import Fraction as Q
from math import comb, factorial as fac

N = 12
SPINS = tuple(range(N, -1, -2))
WEIGHTS = {1: Q(13, 7), 2: Q(7319, 5929), 3: Q(390, 539),
           4: Q(884, 4851), 8: Q(247, 5082)}
BETA = Q(13, 84)
SPIN_WEIGHT = Q(13, 154)


def swap_eigenvalue(s, J):
    """Compressed swap of s qubits on the total-spin J subspace."""
    r = N-J
    c = [Q((-1)**i*comb(r, i), comb(N, i)*comb(N, r-i))
         for i in range(r+1)]
    norm = sum(comb(N, i)*comb(N, r-i)*v*v for i, v in enumerate(c))
    total = Q(0)
    for a in range(s+1):
        for b in range(s+1):
            for i in range(N-s+1):
                j = r-a-b-i
                if 0 <= j <= N-s:
                    total += (comb(s, a)*comb(s, b)*comb(N-s, i)*comb(N-s, j)
                              *c[a+i]*c[b+i])
    return total/norm


def reduction_factor(s, ell):
    return Q(fac(s)**2*fac(N-ell)*fac(N+ell+1),
             fac(N)**2*fac(s-ell)*fac(s+ell+1))


def multipoles(purities):
    weights = []
    for s, purity in enumerate(purities):
        weights.append((purity-sum(reduction_factor(s, ell)*weights[ell]
                                   for ell in range(s)))/reduction_factor(s, s))
    return weights


def verify():
    assert all(w > 0 for w in WEIGHTS.values()) and SPIN_WEIGHT > 0
    for J in SPINS:
        swaps = [swap_eigenvalue(s, J) for s in range(9)]
        weights = multipoles(swaps)
        assert swaps[0] == 1
        assert swaps[6]-BETA == (sum(w*weights[ell] for ell, w in WEIGHTS.items())
                                 + SPIN_WEIGHT*int(J == 2))
    print("Verified the quartic spin-covariant identity on all seven even-spin channels.")

    # c_k = sqrt(binomial(N,k)) z_k / sqrt(norm), using rational z_k.
    z = [Q(0)]*(N+1)
    z[1], z[6], z[11] = Q(1), Q(1, 7), Q(-1)
    norm = sum(comb(N, k)*x*x for k, x in enumerate(z))
    assert [comb(N, k)*z[k]**2/norm for k in (1, 6, 11)] == [Q(7,25),Q(11,25),Q(7,25)]

    def purity(s):
        return sum(comb(s, i)*comb(s, h) *
                   sum(comb(N-s, j)*z[i+j]*z[h+j] for j in range(N-s+1))**2
                   for i in range(s+1) for h in range(s+1))/norm**2

    assert purity(6) == BETA
    assert all(multipoles([purity(s) for s in range(9)])[ell] == 0 for ell in WEIGHTS)
    print("Verified exact attainment of P6 = 13/84 by the icosahedral state.")


if __name__ == '__main__':
    verify()
