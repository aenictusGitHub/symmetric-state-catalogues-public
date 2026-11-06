Exact Dicke coefficient display

This generator reads data/catalogue.json and writes the separate display asset
data/coefficient-latex.js. It maps each unique exactCoefficients expression,
including zero, to SymPy LaTeX. It does not change coefficient values, witnesses,
verification records, or catalogue data. Numerical and implicit algebraic states
without explicit exact coefficients do not acquire exact amplitude formulas.

From the anticoherent-state-catalogue directory, using Python with SymPy:

    python3 research/coefficient-latex/build.py
    python3 research/coefficient-latex/build.py --check

The accepted input grammar is integers, +, -, *, /, parentheses, I, and sqrt
with one argument. Syntax is checked before parsing. Free symbols, floating-point
constants, nonfinite values, and all other operations/functions are rejected.
The standard SymPy LaTeX printer is used without an additional simplify, refine,
expand, or numerical-evaluation step. Output keys retain the original expressions.

Browser API:

    window.AC_COEFFICIENT_LATEX[originalExactExpression]

The mapping is frozen. Render the returned LaTeX with the existing local KaTeX
renderer; the original catalogue remains the source of exact coefficients.
