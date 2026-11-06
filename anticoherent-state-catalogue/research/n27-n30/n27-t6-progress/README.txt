N=27, t=6: rigorous partial results

The unrestricted maximum is now proved to be t_max(27)=5.
The universal six-copy proof is in ../n27-t6/. This folder records
independent auxiliary arguments; it is not that universal certificate.

proof.tex: exact exclusion of proper rotational symmetries of order >=4.
check_symmetry.py: self-contained exact rational checker (400 cases).
symmetry-verification.json: every exact rational witness value.
independent-symmetry-verification.json: independent algebra audit summary.
necessary-condition.tex: conditional Kramers/frame positivity derivations,
  including the rank-four strengthening. These are necessary conditions,
  not a universal exclusion proof.
positive-control.json: exact N49 AC6 quadrature control and floating matrix
  implementation checks, explicitly distinguished from a N27 certificate.
check_positive_control.py: exact rational check for that quadrature state.
search-report.json: scoped numerical searches; no inference from failures.

Run with Python 3 and no third-party dependencies:
  python3 check_symmetry.py
  python3 check_positive_control.py
