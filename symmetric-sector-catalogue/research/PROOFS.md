# Exact two-copy certificates for symmetric sector lengths

These certificates concern normalized pure states in `Sym^n(C²)`, with
`S_k = sum_{wt(P)=k} <P>²` and `S_0=1`. They are new algebraic deductions from the tensor/purity formalism, not additional claims attributed to the supplied draft. All printed numerical optimization is used only to discover a rational certificate. The independent verifier checks every identity using exact fractions.

From the website directory, run `python3 research/certificates/verify_certificates.py`. This uses the standard library only. To reproduce certificate discovery, run `python3 research/certificates/lp_bounds.py` with NumPy, SciPy and SymPy installed. The latter is not required for verification.

## The spectral basis

Let `j=n/2` and `|Ψ>=|ψ>⊗|ψ>`. Decompose the product of its two spin-j registers into total-spin sectors `J=0,...,n`, and set `p_J=<Ψ|Π_J|Ψ>`. Then `p_J≥0`, `sum_J p_J=1`. Interchange of the registers acts as `(-1)^(n-J)`, so only `J=n,n-2,...` occur.

Let `R_t` be the operator swapping the first t qubits between registers, projected onto `Sym^n⊗Sym^n`. Its expectation on `|Ψ>` is the reduced purity `r_t`. It commutes with simultaneous SU(2) rotations, and the total-spin decomposition is multiplicity-free. Hence it has one scalar eigenvalue `h(n,d,t)` on each sector, where `d=n-J`.

For unnormalized Dicke sums `U_a=sum_{|x|=a}|x>`, one has `||U_a||²=C(n,a)` and `J_+ U_a=(n-a+1)U_{a-1}`. A highest-weight vector of total spin `J=n-d` is

`v = sum_{c=0}^d v_c U_c⊗U_{d-c}`,

`v_c = (-1)^c (n-d+c)! (n-c)! / ((n-d)! n!)`.

The coefficients follow directly from `v_{c+1}(n-c)+v_c(n-d+c+1)=0` and `v_0=1`. Taking the `U_0⊗U_d` matrix element of `R_t v` gives

`h(n,d,t) = [sum_{c=0}^d v_c C(t,c) C(n-t,d-c)] / C(n,d)`.

Indeed, to produce a zero-excitation first output register, all c excitations in its input must lie in the swapped block and all d-c excitations of the second input must lie outside that block. The integer count is exactly `C(t,c)C(n-t,d-c)`. This establishes the rational eigenvalue formula without diagonalization or floating-point arithmetic.

By binomial inversion of reduced purities (draft Eq.20), define

`s(n,k,d) = C(n,k) sum_{t=0}^k (-1)^(k+t) 2^t C(k,t) h(n,d,t)`.

Then `S_k=sum_J s(n,k,n-J)p_J`.

The Hilbert-Schmidt multipole length `a_L=sum_M |rho_LM|²` is nonnegative. From draft Eq.22 its eigenvalue in this same basis is

`q(n,L,d) = [(2L+1)C(n,L)/((L+1)C(n+L+1,n))]`
`             × sum_{t=0}^L (-1)^(t+L) C(L,t) C(t+L,t) h(n,d,t)`.

Thus `a_L=sum_J q(n,L,n-J)p_J≥0`. These formulas use the orthonormal multipole convention of the draft, so `a_0=1/(n+1)`.

## How a certificate proves a bound

Each certificate is an exact identity

`B-S_k = sum_i lambda_i g_i` (upper bound), or

`S_k-B = sum_i lambda_i g_i` (lower bound),

where every `lambda_i` is a nonnegative rational number and every `g_i` is known nonnegative on physical pure symmetric states. The verifier expands both sides in the `p_J` basis and checks equality coefficient by coefficient. It also checks that the LP vertex satisfies every basis constraint; such a vertex is not assumed to correspond to a physical state.

The newly sharp bounds below use only `p_J≥0` and `a_L≥0`; no numerical/global-optimality assumptions are involved. The larger interval table also uses:

- `r_t≥1/(t+1)` and `r_t≤1` for `t≤floor(n/2)`;
- `S_2≤C(n,2)` for `n≥3` (published bound cited by draft Table III);
- the published full-sector maximum `S_n≤2^(n-1)` for odd n and `≤2^(n-1)+1` for even n;
- `S_4≤C(n,4)` for `n≥9` (published low-rank bound/maximum threshold in the draft);
- the draft's Theorem 1 lower bound `S_8≥9` at n=9;
- the reduced-state shadow and even-sector inequalities in draft Eq.27, also published in [Serrano-Ensástiga, Giraud and Martin, Physical Review A 113, 012415 (2026)](../references/PRA113_012415(2026).pdf#page=4), PDF p.4, Eqs. (35)–(36). The even-sector inequality is applied to each symmetric reduced state and rescaled by C(n,t).

All dependence on those additional inequalities is visible in each certificate's named terms. An upper/lower bound on a relaxation is not an assertion of attainability.

## Newly sharp certificates

Here `p_J` always means the projection of `|ψ>⊗|ψ>` onto total spin J. All coefficients below multiply nonnegative quantities.

### n=7, maximum S6=49

`49-S6 = (816/7)a2 + (704/7)a4`.

Attained by `(2√2 D7^0 + √14 D7^3 + √14 D7^6)/6`, draft A34.

### n=8, minimum S4=310/21

`S4-310/21 = (2304/7)p0 + (900/7)a1 + 66 a2 + (468/7)a5`.

Attained by `(4 D8^2 + 2 D8^5 + √7 D8^8)/(3√3)`, draft A42.

### n=8, maximum S6=364/3

`364/3-S6 = (594/25)p2 + (7656/25)a2 + (954/5)a5 + (4224/25)a7`.

Attained by `(√30 D8^0 + 2√21 D8^4 + √30 D8^8)/12`, draft A40.

### n=8, maximum S7=592/7

`592/7-S7 = (5328/7)p0 + (1950/7)a1 + (468/7)a5`.

Attained by `(√7 D8^0 + 2 D8^3 + 4 D8^6)/(3√3)`, draft A38.

### n=10, maximum S7=320

`320-S7 = (2288/3)p4 + 2464 p0 + (3344/5)a1 + (1144/5)a3`.

Attained by `(D10^0 + √2 D10^5 + D10^10)/2`, draft A56.

### n=10, maximum S8=405

`405-S8 = (21362/125)p2 + (150238/125)a2 + (233763/325)a5`
`          + (1038906/1625)a7 + (7689/25)a9`.

Attained by `(√5 D10^1 + √6 D10^5 + √5 D10^9)/4`, draft A54.

### n=12, maximum S7=946

`946-S7 = (9310/11)p8 + (11628/11)p6 + (5644/11)p4`
`         + (12896/11)p2 + 7150 p0 + 728 a3`.

Attained by the one-excitation Dicke state `D12^1`. Its sector length is exactly
`C(12,7)[1 + 4·7(2·7-1-12)/12²] = 946`.

### n=12, maximum S10=1606

`1606-S10 = (13267904/27099)p4 + (18264662/3011)a2`
`            + (17261234848/2917659)a7 + (48301344/51187)a8`
`            + (1125641920/972553)a9 + (49898368/57209)a11`.

Attained by the icosahedral representative

`|ψ> = (√7 D12^1 + √11 D12^6 - √7 D12^11)/5`.

### n=10,11,12, maximum S5=C(n,5)

The spectral certificates contain only nonnegative pJ terms. The exact coefficients are in `sharp-certificates.json`. The coherent state `D_n^0` lies entirely in J=n, where the eigenvalue is `C(n,5)`, and therefore attains the bound. In particular, the ten-qubit bound plus convexity and symmetric reduction extends it to every n≥10.

## Remaining gaps

`lp-bounds.json` gives 154 verified lower/upper bounds for every sector at n=2,...,12. Only a bound with a matching physical witness is a solved extremum. For example, this relaxation gives `20/3≤S6` at n=7, not the proposed minimum 7; at n=9 it gives `S6≥36`, not the candidate 420/11. It gives `S6≤8348/9` at n=12, leaving a gap above the coherent value 924. These cases are not proved by these certificates.
