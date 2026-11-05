# No pure 4-anticoherent state of eleven symmetric qubits

**Theorem.** There is no nonzero pure state in Sym^11(C²) that is anticoherent to order 4. The result covers arbitrary complex Dicke coefficients; no discrete symmetry, sparse support, or reality assumption is imposed.

## A scalar polynomial identity

Let |ψ⟩ be an arbitrary, not necessarily normalized spin-11/2 state. Choose irreducible tensor operators T_{lq}, with l=0,…,11 and q=−l,…,l. Define the homogeneous cubic invariant in their expectation values

I_{lmk}(ψ) = Σ_{q+r+s=0} ( l m k ; q r s ) ⟨ψ|T_{lq}|ψ⟩ ⟨ψ|T_{mr}|ψ⟩ ⟨ψ|T_{ks}|ψ⟩,

where the parenthesis denotes a Wigner 3j symbol. Only triples obeying the triangle inequalities and with l+m+k even occur below. Normalize each invariant at the coherent state |D_11^(0)⟩:

W_{lmk}(ψ) = I_{lmk}(ψ) / I_{lmk}(D_11^(0)).

All these denominators are nonzero: each coherent tensor moment with q=0 and l≤11 is nonzero, and the 3j symbol with q=r=s=0 is nonzero for the listed even triangle triples. Scaling every rank-l tensor family by any nonzero rank-dependent constant does not change W.

The exact identity is

Σ_{lmk} b_{lmk} W_{lmk}(ψ) = ⟨ψ|ψ⟩³.

The coefficient table follows. The final column is the equivalent rational coefficient α used by the verifier's tensor convention, defined below.

| (l,m,k) | b for coherent-normalized W | α for verifier invariant f |
|---|---:|---:|
| (0,1,1) | 9 | -18/121 |
| (0,2,2) | 56625/2303 | 13590/278663 |
| (0,3,3) | 4779/329 | -2124/199045 |
| (0,4,4) | 3204/2303 | 178/199045 |
| (1,1,2) | 180/7 | 216/9317 |
| (1,2,3) | 8910/329 | -216/39809 |
| (1,3,4) | -4272/329 | -1424/938355 |
| (1,4,5) | -1560/329 | 312/437899 |
| (2,2,2) | 27500/16121 | -120/1950641 |
| (2,2,4) | -295650/16121 | -3942/3065293 |
| (2,3,3) | -4320/329 | -192/2189495 |
| (2,3,5) | -103500/3619 | 8280/4816889 |
| (2,4,4) | -1091700/177331 | 3639/168591115 |
| (2,4,6) | -101250/25333 | -1350/3065293 |
| (2,5,5) | -900/4277 | -36/39848809 |
| (2,5,7) | 4050/47047 | -108/3065293 |
| (3,3,4) | -3618/329 | 268/4691775 |
| (3,3,6) | -23400/3619 | -208/437899 |
| (3,4,5) | -20520/47047 | -114/62619557 |
| (3,4,7) | -10050/6721 | 268/938355 |
| (3,5,6) | 450/517 | -20/3065293 |
| (4,4,4) | 1899774/2305303 | 3909/3130977850 |
| (4,4,6) | 26280/25333 | -146/15326465 |
| (4,4,8) | -360/6721 | -12/312785 |
| (4,5,5) | 6750/4277 | -75/39848809 |
| (4,5,7) | 9600/47047 | 32/6568485 |

Every triple in the table contains a rank between 1 and 4. Four-anticoherence means every tensor moment of ranks 1,2,3,4 vanishes. Thus every summand on the left vanishes, while the right is strictly positive for a nonzero state. This is the contradiction.

## Exact rational convention used by the verifier

Use the weighted binary-form basis e_r=√C(11,r)|D_11^(r)⟩ and write |ψ⟩=Σ_r z_r e_r. Its Gram matrix is diagonal, G_rr=C(11,r). We explicitly choose the spin-axis convention J_z|D_11^(r)⟩=(r−11/2)|D_11^(r)⟩. This reverses the source convention m=11/2−r; anticoherence and the scalar invariant are rotation-invariant, so this convention change does not affect the conclusion. With this choice, angular momentum acts as

J_+ e_r=(r+1)e_(r+1),     J_- e_r=(12−r)e_(r−1).

Define F_{lq}=(ad J_-)^(l−q)(J_+)^l / [(l−q)! l!], where ad A(B)=AB−BA. These matrices have integer entries in this basis. For a fixed rank triple, define

h_{lmk}^{qrs} = ( l m k ; q r s ) / √[C(2l,l−q) C(2m,m−r) C(2k,k−s)].

Let h_seed be its first nonzero value in ascending lexicographic (q,r) order, with s=−q−r. Set g=h/h_seed. This denominator is nonzero by construction; each resulting g is rational. The verifier invariant is

f_{lmk}(ψ)=Σ_{q+r+s=0} g_{lmk}^{qrs} ⟨F_{lq}⟩ ⟨F_{mr}⟩ ⟨F_{ks}⟩.

The binomial factors are the precise conversion from the commutator-normalized F tensors to a conventional spherical-tensor family: F_{lq} is proportional, by a constant depending only on l, to √C(2l,l−q) T_{lq}. Thus f is a nonzero scalar multiple of I. Dividing by the coherent value gives exactly the same W.

In this convention, f(D_11^(0))=g^{000} C(11,l)C(11,m)C(11,k), and b=α f(D_11^(0)). Both relations and Σb=1 are verified exactly. The α coefficients certify Σαf=⟨ψ|ψ⟩³.

## Reproduction and scope

Run `python3 verify.py` in this directory (Python 3 and SymPy are the only dependencies). The script constructs every F matrix from actual commutators, expands the full polynomial in independent formal variables z_r and conjugate(z_r), and checks every coefficient against (Σ_r C(11,r)|z_r|²)³. It uses no floating-point arithmetic, numerical fit, highest-weight decomposition, or optimization output as evidence.

The independent expansion checks 974,308 contributions over 5,664 monomial slots, with zero nonzero residual coefficients. The exact identity was also separately obtained by rational highest-weight block algebra on Sym³(C¹²), but that decomposition is not needed to trust this verifier.

This is an independently derived certificate for this catalogue. It resolves the previously unresolved N=11, t=4 catalogue cell. An exact third-order witness is already known, so the maximal anticoherence order for eleven qubits is exactly three. The standard multipole and reduced-state characterizations of t-anticoherence are given in PRA 92, 052333 (2015), Eqs. (10), (16)–(17), PDF pages 3–4.
