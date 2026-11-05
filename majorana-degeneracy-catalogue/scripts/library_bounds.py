"""Rigorous partition bounds; numerical reported AC orders are never inputs.

Proof ingredients already verified in research/SMALL_N_PROOFS.md and
research/partition-proofs.md. Suitable for source-observed partitions up to
N=100 without enumerating the other partitions.
"""
from math import isqrt

CORE_GLOBAL_UPPER = {2:1, 3:1, 4:2, 5:1, 6:3, 7:2, 8:3,
                     9:2, 10:3, 11:3, 12:5}


def rigorous_order_upper(n, partition):
    """Upper bound for the maximal AC order in this exact full partition."""
    lam = tuple(sorted(partition, reverse=True))
    if not isinstance(n, int) or n < 2 or not lam:
        raise ValueError('Require integer n>=2 and a nonempty partition')
    if any(not isinstance(x, int) or x <= 0 for x in lam) or sum(lam) != n:
        raise ValueError('Partition entries must be positive integers summing to n')
    g = lam[0]
    balanced = n % 2 == 0 and lam == (n // 2, n // 2)
    # Complete first-order criterion, including the non-polystable boundary.
    if 2*g > n or (2*g == n and not balanced):
        return 0
    # Schmidt rank and the established tensor Gram bound for q=floor(t/2).
    upper = min(n // 2, 2*isqrt(n + 1) - 1)
    # At n=2t+1 and t>=2, a complementary rank-one kernel is impossible.
    if n >= 5 and n % 2:
        upper = min(upper, (n - 3) // 2)
    # Previously proved small-N global maxima; not numerical library maxima.
    if n in CORE_GLOBAL_UPPER:
        upper = min(upper, CORE_GLOBAL_UPPER[n])
    # Existing degree-four certificate at N12 rules out every double point.
    if n == 12 and g >= 2:
        upper = min(upper, 3)
    # Degree-two support moment: E[(K-g)(n-K)] = n(n-1-3g)/6.
    if 3*g > n - 1:
        upper = min(upper, 1)
    # Degree-three moment matrix: (n/2-g)^2 >= n(n+2)/12.
    if 3*(n - 2*g)**2 < n*(n + 2):
        upper = min(upper, 2)
    # Equality in the second-order bound fixes the entire partition.
    if n == 3*g + 1 and lam != (g,) + (1,)*(n - g):
        upper = min(upper, 1)
    return upper


def theorem_lower(n, partition):
    """Guaranteed lower bound without importing any numerical witness claim."""
    upper = rigorous_order_upper(n, partition)
    return int(upper > 0)  # finite-filter theorem or balanced Dicke witness


PROOF_IDS = {
    'firstOrder': 'partition-first-order-exclusion',
    'firstOrderExistence': 'partition-filtering-existence',
    'schmidt': 'ac-schmidt-rank',
    'gram': 'ac-tensor-gram-even',
    'oddComplement': 'ac-odd-complement-rank',
    'secondMoment': 'order-two',
    'thirdMoment': 'order-three',
    'secondOrderSaturation': 'partition-second-order-saturation',
    'n12RepeatedPoint': 'n12-higher',
}


if __name__ == '__main__':
    for n in (12, 13, 20, 24, 42, 100):
        print(n, '(1^n)', theorem_lower(n, [1]*n), rigorous_order_upper(n, [1]*n))
