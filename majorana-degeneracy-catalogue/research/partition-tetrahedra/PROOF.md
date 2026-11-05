# Two exact twelve-qubit partition witnesses

Let T_+(z)=z^4+2 i sqrt(3) z^2+1 and T_-(z)=z^4-2 i sqrt(3) z^2+1. Each has four simple roots; their greatest common divisor is1. These claims are checked exactly in the verifier.

Use the Majorana convention P(1,z)=sum_r (-1)^r sqrt(binomial(12,r)) c_r z^r, up to an overall nonzero scalar. For a monic polynomial P, take c_r proportional to (-1)^r [z^r]P/sqrt(binomial(12,r)) and normalize exactly.

For P=T_+^3 the complete partition is (3,3,3,3). For P=T_+^2 T_- it is (2,2,2,2,1,1,1,1). There is no point at infinity in either case because both affine polynomials have degree12.

The exact normalized Dicke coefficients are provided in states.json. Direct evaluation gives rho_2=I_3/3 for both states, establishing order2 in both full-partition cells. Their squared third-order reduced-matrix distances from I_4/4 are441/14884 and1/484, respectively, so the displayed representatives have exact order2.

The partition(3,3,3,3) has maximum order2 once combined with the independently established multiplicity-three exclusion of order3 at n12. The second representative supplies a lower bound2 for partition(2,2,2,2,1,1,1,1); its failure at order3 does not exclude other order3 representatives in the same partition.

Run `python3 verify.py` from any directory. It requires SymPy, reads the neighboring states.json, and uses exact arithmetic throughout.
