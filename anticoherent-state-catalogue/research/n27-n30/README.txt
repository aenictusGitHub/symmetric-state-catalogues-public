EXACT PURE AC ATTAINMENT: N=27--30

radical-witnesses/states.json contains explicit normalized radical Dicke
vectors for N27 order5, N29 order5 and N30 order7. Run
  python3 radical-witnesses/verify_witnesses.py
with SymPy installed. The verifier reads the saved exact vectors and computes
direct partial traces independently of build_witnesses.py. It also proves
that these states fail at the next order and verifies Majorana multiplicities.

n28-t6/certificate.json defines a unique real algebraic root near the refined
Liège acs.28.6 state. It is an exact implicit definition, not a list of radical
Dicke components. Run
  python3 n28-t6/verify_independently.py
with standard Python3. It reconstructs all rational equations, Jacobian and
Banach bounds without numerical dependencies. A rational interval also proves
this root is not seventh-order AC.

build_states.py combines the four records for catalogue import and generates
numerical sphere coordinates from the radical states. It requires SymPy,
mpmath and NumPy. When this folder is installed as research/n27-n30, its
source defaults to data/catalogue.json; alternatively pass --source PATH.
Only N28 geometry is copied from the nearby numerical library representative.
The exact N29 construction has two double poles and25 simple finite roots;
its partition is (2,2,1^25), unlike an all-simple numerical source.

The constructions establish lower bounds. Global upper exclusions are
separate artifacts supplied by the parent workflow. In particular, N30
order8 has been independently excluded, so its exact order7 vector attains
the global maximum. The combined-state builder reads proofs.json when available and marks a
global maximum only when its exact attaining order meets a proved upper bound.

proof.tex in each subfolder states the construction and proof scope.

Updated universal bounds: t_max(27) and t_max(29) are both in[5,6].
N28 has proved maximum6, and N30 has proved maximum7.
