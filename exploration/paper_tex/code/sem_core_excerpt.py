# -*- coding: utf-8 -*-
"""Excerpt 1: graded spectral-element assembly (from sem_core.py).

Even-Legendre modal basis  phi_j(sigma) = P_{2j}(2*sigma-1)  on sigma in [0,1],
Gauss-Legendre quadrature with the cylindrical measure  sigma dsigma.
The mass matrix is EXACTLY diagonal:  M_jj = 1 / (2*(4j+1)).
"""
import numpy as np
from numpy.polynomial import polynomial as P
from scipy.special import roots_legendre


def cheb_lobatto(n):
    """n+1 Chebyshev-Lobatto nodes on [-1,1] (ascending)."""
    return -np.cos(np.pi * np.arange(n + 1) / n)


class SEM1D:
    def __init__(self, breaks=(0.0, 0.6, 0.9, 1.0), n=18, nq=None):
        self.n, self.ne = int(n), len(breaks) - 1
        self.nq = int(nq) if nq else self.n + 3
        self.elems, gid = [], 0
        for e in range(self.ne):
            a, b = breaks[e], breaks[e + 1]
            he = b - a
            xi = cheb_lobatto(self.n)                     # element nodes
            # Lagrange polynomials of the element nodes (coefficients in xi)
            Lc = []
            for k in range(self.n + 1):
                roots_m = np.delete(xi, k)
                Lc.append(P.polyfromroots(roots_m)
                          / np.prod(xi[k] - roots_m))
            # Gauss-Legendre quadrature nodes mapped to [a,b]
            xq, wq = roots_legendre(self.nq)
            sig_q = 0.5 * (a + b) + 0.5 * he * xq
            w_sig = 0.5 * he * wq * sig_q                 # measure sigma dsigma
            L = np.array([P.polyval(xq, c) for c in Lc])                  # (n+1,nq)
            DL = np.array([P.polyval(xq, P.polyder(c)) * (2.0 / he)
                           for c in Lc])                                  # d/dsigma
            g = (np.arange(self.n + 1) if e == 0 else
                 np.concatenate([[gid - 1], np.arange(gid, gid + self.n)]))
            gid = gid + self.n if e else self.n + 1
            self.elems.append(dict(a=a, b=b, he=he, xi=xi, sig_q=sig_q,
                                   w=w_sig, L=L, DL=DL, g=g))
        self.ndof = gid
        self.sig_nodes = np.zeros(self.ndof)
        for el in self.elems:
            self.sig_nodes[el['g']] = 0.5 * (el['a'] + el['b']) \
                + 0.5 * el['he'] * el['xi']
        self.surf_dof = int(np.argmax(self.sig_nodes))     # node at sigma = 1

    def gather(self, u):
        """Nodal values -> values at all quadrature points."""
        return [el['L'].T @ u[el['g']] for el in self.elems]

    def assemble(self, coefs, kind='L'):
        """Global matrix  sum_q w_q * coef_q * (L or DL)[:,q] (.)^T."""
        out = np.zeros((self.ndof, self.ndof))
        for el, cf in zip(self.elems, coefs):
            B = el['L'] if kind == 'L' else el['DL']
            out[np.ix_(el['g'], el['g'])] += (B * (el['w'] * cf)) @ B.T
        return out
