# An exact algebraic eleven-qubit third-order AC state with a double Majorana point

The state is specified exactly as an isolated solution of rational quadratic equations in a rational box. Decimal amplitudes are provided for display, but the displayed decimals themselves are not asserted to satisfy the equations exactly.

## Exact definition of the state

Write

  c_r = sqrt(binomial(11,r)) (x_r + i y_r), 0 <= r <= 9,
  c_10 = c_11 = 0.

Order the twenty real coordinates as x_0,...,x_9,y_0,...,y_9. Fix coordinates 0,10,12,19 to the rational values given in witness.json. The other sixteen coordinates lie in the closed infinity-norm box of radius 10^(-30) centered at the exact rational entries stored in that file.

For a,b=0,...,3, define

  G_ab = sum_(q=0)^8 binomial(8,q) z_(a+q) conjugate(z_(b+q)),
  z_r=x_r+i y_r for r<=9, and z_10=z_11=0.

The exact defining equations are

  G_aa=1/[4 binomial(3,a)],
  Re G_ab=Im G_ab=0 for a<b.

These are sixteen real quadratic equations F=0 with rational coefficients. They imply rho_3=I_4/4 because rho_3[a,b]=sqrt(binomial(3,a)binomial(3,b)) G_ab. They therefore imply exact normalization and third-order anticoherence.

## Existence, uniqueness and algebraicity

Let x0 be the rational center, J0=DF(x0), R=J0^(-1), and define T(x)=x-RF(x). All quantities in the following verification are calculated exactly over the rational numbers. Let

  e = ||R F(x0)||_infinity,
  M = max_a sum_i |R_ai| sum_(j,k) |d_j d_k F_i|.

The Hessians are constant because the equations are quadratic. The verifier reconstructs every equation and checks

  e < 10^(-70),
  M < 12,000,000.

For r=10^(-30), Taylor's formula gives

  ||T(x)-x0||_infinity <= e + M r^2 < r,
  ||DT(x)||_infinity <= M r < 1.2 * 10^(-23)

throughout the box. The Banach fixed-point theorem therefore gives a unique real root in the sixteen-variable slice of the box with the four specified coordinates fixed. The derivative stays invertible there, so this root is isolated and algebraic. This proves an exact state exists; it does not infer exactness from a small numerical residual.

The exact rational center was obtained by refining the rounded example in PRA96_032304(2017), PDF p.11, Appendix A.1, Eq. (A1), after reversing its spin-basis order to our Dicke excitation convention. The existence certificate and algebraic specification are independently supplied here, rather than attributed to the printed decimals.

## Full Majorana multiplicity partition

Use the homogeneous binary polynomial

  P(X,Y)=sum_(r=0)^11 (-1)^r sqrt(binomial(11,r)) c_r X^(11-r)Y^r.

Its affine form f(z)=P(1,z) has coefficients (-1)^r binomial(11,r)(x_r+i y_r). Since c_10=c_11=0 and c_9 is nonzero everywhere in the certified box, f has degree nine. The projective point [0:1] is therefore a root of multiplicity exactly two.

All nine finite roots are simple. This is checked without numerical root clustering. Build the 17-by-17 Sylvester matrix S of f and f'. A rational complex matrix R_S, stored in majorana-certificate.json, satisfies throughout the entire state box

  ||I-R_S S||_infinity < 5.77 * 10^(-22) < 1.

The verifier computes a rigorous rational bound using |Re z|+|Im z| for each complex entry. Hence S is invertible, f and f' are coprime, and the finite roots are simple. The exact Majorana partition is (2,1,1,1,1,1,1,1,1,1), with ten distinct points and maximal multiplicity two.

## Reproduction

Run `python3 verify.py` and `python3 verify_majorana.py`. Both resolve neighboring data paths independently of the current directory. The first requires SymPy and takes about one second; the second uses only the Python standard library. Neither uses a numerical solver or floating-point arithmetic in its proof checks.

Required files: witness.json, verify.py, majorana-certificate.json, verify_majorana.py. The recorded verification JSON files are outputs, not trusted inputs. build_witness.py and build_majorana.py are optional discovery scripts and are unnecessary for verification.
