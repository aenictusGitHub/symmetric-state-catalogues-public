# No pure 4-anticoherent state of ten symmetric qubits

**Theorem.** There is no nonzero pure state in Sym^10(C²) that is anticoherent to order 4. The result covers arbitrary complex Dicke coefficients; no discrete symmetry, sparse support, or reality assumption is imposed.

## A scalar polynomial identity

Let |ψ⟩ be an arbitrary, not necessarily normalized spin-5 state. Choose irreducible tensor operators T_{lq}, with l=0,…,10 and q=−l,…,l. Define the homogeneous cubic invariant in their expectation values

I_{lmk}(ψ) = Σ_{q+r+s=0} ( l m k ; q r s ) ⟨ψ|T_{lq}|ψ⟩ ⟨ψ|T_{mr}|ψ⟩ ⟨ψ|T_{ks}|ψ⟩,

where the parenthesis denotes a Wigner 3j symbol. Only triples obeying the triangle inequalities and with l+m+k even occur below. Normalize each invariant at the coherent state |D_10^(0)⟩:

W_{lmk}(ψ) = I_{lmk}(ψ) / I_{lmk}(D_10^(0)).

All these denominators are nonzero: each coherent tensor moment with q=0 and l≤10 is nonzero, and the 3j symbol with q=r=s=0 is nonzero for the listed even triangle triples. Scaling every rank-l tensor family by any nonzero rank-dependent constant does not change W.

The exact identity is

Σ_{lmk} b_{lmk} W_{lmk}(ψ) = ⟨ψ|ψ⟩³.

The coefficient table follows. The final column is the equivalent rational coefficient α used by the verifier's tensor convention, defined below.

| (l,m,k) | b for coherent-normalized W | α for verifier invariant f |
|---|---:|---:|
| (0,1,1) | 115/4 | -23/40 |
| (0,2,2) | 225/4 | 1/6 |
| (0,3,3) | 135/4 | -3/64 |
| (0,4,4) | 37/4 | 37/2520 |
| (1,1,2) | 60 | 2/25 |
| (1,2,3) | 135/2 | -1/40 |
| (1,3,4) | 5/6 | 1/4320 |
| (1,4,5) | -275/12 | 11/1008 |
| (2,2,4) | -25/2 | -1/486 |
| (2,3,3) | -135/4 | -1/1920 |
| (2,3,5) | -125/2 | 5/432 |
| (2,4,4) | -3175/132 | 127/498960 |
| (2,4,6) | -875/44 | -1/108 |
| (3,3,4) | -1495/44 | 299/570240 |
| (3,3,6) | -625/44 | -5/1152 |
| (3,4,5) | -425/26 | -85/314496 |
| (3,4,7) | -6125/858 | 7/864 |
| (4,4,4) | -279/52 | -31/982800 |
| (4,4,6) | -25/12 | 11/105840 |
| (4,4,8) | -175/286 | -1/252 |

Every triple in the table contains a rank between 1 and 4 (in fact, at least two such ranks). Four-anticoherence means every tensor moment of ranks 1,2,3,4 vanishes. Thus every summand on the left vanishes, while the right is strictly positive for a nonzero state. This is the contradiction.

## Exact rational convention used by the verifier

Use the weighted binary-form basis e_r=√C(10,r)|D_10^(r)⟩ and write |ψ⟩=Σ_r z_r e_r. Its Gram matrix is diagonal, G_rr=C(10,r). For this calculation choose the spin-axis convention J_z|D_10^(r)⟩=(r−5)|D_10^(r)⟩. This reverses the catalogue's default quantization axis; anticoherence and the scalar identity are unchanged. In this convention angular momentum acts as

J_+ e_r=(r+1)e_(r+1),     J_- e_r=(11−r)e_(r−1).

Define F_{lq}=(ad J_-)^(l−q)(J_+)^l / [(l−q)! l!], where ad A(B)=AB−BA. These matrices have integer entries in this basis. For a fixed rank triple, define

h_{lmk}^{qrs} = ( l m k ; q r s ) / √[C(2l,l−q) C(2m,m−r) C(2k,k−s)].

Let h_seed be its first nonzero value in ascending lexicographic (q,r) order, with s=−q−r. Set g=h/h_seed. This denominator is nonzero by construction; each resulting g is rational. The verifier invariant is

f_{lmk}(ψ)=Σ_{q+r+s=0} g_{lmk}^{qrs} ⟨F_{lq}⟩ ⟨F_{mr}⟩ ⟨F_{ks}⟩.

The binomial factors are the precise conversion from the commutator-normalized F tensors to a conventional spherical-tensor family: F_{lq} is proportional, by a constant depending only on l, to √C(2l,l−q) T_{lq}. Thus f is a nonzero scalar multiple of I. Dividing by the coherent value gives exactly the same W.

In this convention, f(D_10^(0))=g^{000} C(10,l)C(10,m)C(10,k), and b=α f(D_10^(0)). Both relations and Σb=1 are verified exactly. The α coefficients certify Σαf=⟨ψ|ψ⟩³.

## Reproduction and scope

Run `python3 verify.py` in this directory (Python 3 and SymPy are the only dependencies). The script constructs every F matrix from actual commutators, expands the full polynomial in independent formal variables z_r and conjugate(z_r), and checks every coefficient against (Σ_r C(10,r)|z_r|²)³. It uses no floating-point arithmetic, numerical fit, highest-weight decomposition, or optimization output as evidence.

The independent expansion checks 430,712 contributions over 3,790 monomial slots, with zero nonzero residual coefficients. The exact identity was also separately obtained by rational highest-weight block algebra on Sym³(C¹¹), but that decomposition is not needed to trust this verifier.

This is an independently derived certificate for this catalogue. Prior published numerical/nonexistence reports motivated the question but are not used as proof. The standard equivalence between t-anticoherence and vanishing multipoles of ranks 1,…,t is stated in the supplied manuscript 2605.29436.pdf, Definition 1, Eqs. (3)–(4), PDF page 2; the reduced-state characterization is also given in PRA 92, 052333 (2015), Eqs. (10), (16)–(17), PDF pages 3–4.
