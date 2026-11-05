Twenty-five qubits: maximum pure AC order five

Unpack the complete bundle and run:
  python3 verify_independently.py

Requires only Python 3 and its standard library. The complete independent
verification typically takes a few minutes. No optimization package or
saved discovery cache is used.

The verifier proves universal nonexistence of pure sixth-order AC states
for N=25, by exact rational arithmetic. Exact fifth-order attainment is
verified separately in ../n25-t5/. See proof.tex for the full argument.

Optional physical controls:
  python3 positive_control.py
  python3 verify_recoupling_control.py
