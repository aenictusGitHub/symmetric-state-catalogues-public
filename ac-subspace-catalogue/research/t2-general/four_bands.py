#!/usr/bin/env python3
"""Universal four-band AC2 construction and exact all-N polynomial certificate.

K(N,2) >= floor(N/8)-3 for N >= 40. This is not an optimality claim.
The coefficient certificate is symbolic for all N, not a finite N scan.
"""
from pathlib import Path
import json
import sympy as S
import numpy as np


def symbolic_weights(x, v):
    return [v/4+x*x/4-S.Rational(35,36)*x+S.Rational(5,6),
            -v/4-x*x/4+S.Rational(5,12)*x+S.Rational(1,2),
            -v/4-x*x/4+S.Rational(13,12)*x-S.Rational(1,2),
            v/4+x*x/4-S.Rational(19,36)*x+S.Rational(1,6)]


def universal_certificate():
    q, i, z, x, v = S.symbols('q i z x v', nonnegative=True)
    weights = symbolic_weights(x, v)
    assert S.expand(sum(weights)-1) == 0
    assert S.expand(sum(weights[l]*(x-l) for l in range(4))) == 0
    assert S.expand(sum(weights[l]*(x-l)**2 for l in range(4))-v) == 0
    records = []
    for remainder in range(8):
        n = 8*q+remainder
        start = 4+remainder//2
        j, h = n/2, 2*q
        xx = (j-start-2*i)/h
        vv = j*(j+1)/(3*h*h)
        # Clearing this common positive denominator produces integer polynomials.
        W = [S.cancel(576*q*q*w) for w in symbolic_weights(xx, vv)]
        assert all(S.denom(w) == 1 for w in W)
        weight_polys = [S.Poly(S.expand(w.subs(q,i+z+4)),i,z) for w in W]
        assert all(c >= 0 for w in weight_polys for c in w.coeffs())
        assert all(w.coeff_monomial(1)>0 for w in weight_polys)
        E = []
        for ell in range(4):
            rr = start+2*i+2*q*ell
            e = W[ell]*W[ell].subs(i,i+1)*(rr+1)*(rr+2)*(n-rr)*(n-rr-1)
            E.append(S.Poly(S.expand(e.subs(q,i+z+5)),i,z))
        assert all(c >= 0 for e in E for c in e.coeffs())
        F = sum(E)**2-2*sum(e*e for e in E)
        A = S.Poly(7421703487488*(i**14*z**2+i**2*z**14),i,z)
        positive_F = F+A
        positive_product = 64*E[0]*E[1]*E[2]*E[3]-A*A
        assert all(c >= 0 for c in positive_F.coeffs())
        assert all(c >= 0 for c in positive_product.coeffs())
        assert positive_F.coeff_monomial(1)>0
        assert positive_product.coeff_monomial(1)>0
        records.append({
            'remainder_mod_8': remainder,
            'all_row_weight_substitution': 'q=i+z+4',
            'edge_substitution': 'q=i+z+5',
            'weight_terms': [len(w.terms()) for w in weight_polys],
            'weight_positive_constants': [str(w.coeff_monomial(1)) for w in weight_polys],
            'F_plus_A_terms': len(positive_F.terms()),
            'F_plus_A_positive_constant': str(positive_F.coeff_monomial(1)),
            '64_product_E_minus_A_squared_positive_constant': str(positive_product.coeff_monomial(1)),
            '64_product_E_minus_A_squared_terms': len(positive_product.terms()),
            'all_coefficients_nonnegative': True,
        })
    return {'theorem': 'K(N,2) >= floor(N/8)-3 for every integer N>=40',
            'arithmetic': 'symbolic exact integer and rational',
            'certificate_type': 'All-N coefficient positivity in two nonnegative variables',
            'edge_scale': 'E_ell=(576 q^2)^2 L_ell^2',
            'A': '7421703487488*(i^14*z^2+i^2*z^14)',
            'moment_identities_verified': True,
            'remainders': records}


def numerical_basis(n):
    """Numerically render the exact radical prescription; no optimization."""
    if n < 40:
        raise ValueError('This family is stated for N >= 40.')
    q, remainder = divmod(n,8)
    k, h, start = q-3, 2*q, 4+remainder//2
    j = n/2
    v = j*(j+1)/(3*h*h)
    supports = start+2*np.arange(k)[:,None]+h*np.arange(4)[None,:]
    xx = (j-start-2*np.arange(k))/h
    W = np.asarray([[float(w) for w in symbolic_weights(x,v)] for x in xx])
    assert np.min(W)>0
    phases = np.ones((k,4), dtype=complex)
    margins = []
    for i in range(k-1):
        rr = supports[i].astype(float)
        lengths = np.sqrt(W[i]*W[i+1]*(rr+1)*(rr+2)*(n-rr)*(n-rr-1))
        a,b,c,d = lengths
        lower, upper = max(abs(a-b),abs(c-d)), min(a+b,c+d)
        assert upper>=lower-1e-10
        diagonal = (lower+upper)/2
        def pair(u,v):
            cu = (diagonal**2+u*u-v*v)/(2*diagonal*u)
            cv = (diagonal**2+v*v-u*u)/(2*diagonal*v)
            return (cu+1j*np.sqrt(max(0,1-cu*cu)),
                    cv-1j*np.sqrt(max(0,1-cv*cv)))
        z0,z1 = pair(a,b)
        z2,z3 = pair(c,d)
        increments = np.array([z0,z1,-z2,-z3])
        assert abs(np.dot(lengths,increments))<1e-7
        phases[i+1] = phases[i]*increments
        margins.append(float(sum(lengths)-2*max(lengths)))
    basis = np.zeros((n+1,k),dtype=complex)
    for i in range(k):
        basis[supports[i],i] = np.sqrt(W[i])*phases[i]
    return basis, supports, W, margins


def direct_check(n):
    B, supports, weights, margins = numerical_basis(n)
    j = n/2
    Jz = np.diag(j-np.arange(n+1))
    Jp = np.diag(np.sqrt(np.arange(1,n+1)*(n+1-np.arange(1,n+1))),1)
    v = j*(j+1)/3
    blocks = [B.conj().T@B-np.eye(B.shape[1]), B.conj().T@Jz@B,
              B.conj().T@(Jz@Jz)@B-v*np.eye(B.shape[1]),
              B.conj().T@Jp@B, B.conj().T@(Jp@Jp)@B,
              B.conj().T@(Jz@Jp+Jp@Jz)@B]
    residual = max(float(np.max(np.abs(block))) for block in blocks)
    assert residual<1e-8
    return {'n':n,'k':B.shape[1],'max_compression_residual':residual,
            'minimum_quadrilateral_margin':min(margins)}


if __name__ == '__main__':
    directory=Path(__file__).resolve().parent
    certificate=universal_certificate()
    certificate['numerical_rendering_checks']=[direct_check(n) for n in range(40,81)]
    (directory/'four-band-certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    print('Universal exact coefficient certificate passed for all 8 remainders.')
    print('Direct spin compressions checked for N=40,...,80.')
