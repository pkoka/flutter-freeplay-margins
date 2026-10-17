"""
Self-contained continuous-time H-infinity synthesis (Doyle, Glover, Khargonekar,
Francis 1989), two Riccati equations, no dependency beyond scipy.

Generalized plant, simplified regular case:
    xdot = A x + B1 w + B2 u
    z    = C1 x + D12 u          (D11 = 0)
    y    = C2 x + D21 w          (D22 = 0)
with the standard normalisations  D12' D12 = I,  D12' C1 = 0,
                                  D21 D21' = I,  B1 D21' = 0.

For a given gamma the central controller exists iff
    X >= 0 solves   A'X + XA + X(g^-2 B1 B1' - B2 B2')X + C1'C1 = 0
    Y >= 0 solves   AY + YA' + Y(g^-2 C1'C1 - C2'C2)Y + B1 B1' = 0
    rho(XY) < gamma^2
"""
import numpy as np
from scipy.linalg import solve_continuous_are, inv, eigvals


def _are(A, G, Q):
    """Solve A'X + XA - X G X + Q = 0 with G possibly indefinite, via the
    Hamiltonian. scipy's solver assumes G >= 0, which fails here because the
    gamma term makes G indefinite, so take the stable invariant subspace."""
    n = A.shape[0]
    H = np.block([[A, -G], [-Q, -A.T]])
    w, V = np.linalg.eig(H)
    idx = np.argsort(w.real)[:n]              # stable subspace
    if np.max(w.real[idx]) >= -1e-9:
        return None
    U = V[:, idx]
    U1, U2 = U[:n], U[n:]
    if np.linalg.matrix_rank(U1, tol=1e-9) < n:
        return None
    X = np.real(U2 @ np.linalg.inv(U1))
    return 0.5*(X + X.T)


def hinf_central(A, B1, B2, C1, C2, D12, D21, gam):
    """Return (Ac, Bc, Cc, Dc) of the central controller, or None."""
    n = A.shape[0]
    g2 = gam**-2
    X = _are(A, B2 @ B2.T - g2*(B1 @ B1.T), C1.T @ C1)
    if X is None or np.min(np.linalg.eigvalsh(X)) < -1e-7:
        return None
    Y = _are(A.T, C2.T @ C2 - g2*(C1.T @ C1), B1 @ B1.T)
    if Y is None or np.min(np.linalg.eigvalsh(Y)) < -1e-7:
        return None
    if max(abs(eigvals(X @ Y))) >= gam**2:
        return None
    F = -B2.T @ X
    L = -Y @ C2.T
    Z = inv(np.eye(n) - g2*(Y @ X))
    Ac = A + g2*(B1 @ B1.T @ X) + B2 @ F + Z @ L @ C2
    return Ac, -Z @ L, F, np.zeros((F.shape[0], C2.shape[0]))


def hinf_syn(A, B1, B2, C1, C2, D12, D21, lo=1e-2, hi=1e6, tol=1e-3):
    """Gamma bisection. Returns (K, gamma) or (None, inf)."""
    best = None
    if hinf_central(A, B1, B2, C1, C2, D12, D21, hi) is None:
        return None, np.inf
    for _ in range(60):
        mid = np.sqrt(lo*hi)
        K = hinf_central(A, B1, B2, C1, C2, D12, D21, mid)
        if K is None:
            lo = mid
        else:
            best, hi = (K, mid), mid
        if hi/lo < 1 + tol:
            break
    return (best[0], best[1]) if best else (None, np.inf)


def lft_closed_loop(A, B, C, K):
    """Plant (A,B,C,D=0) with controller K=(Ac,Bc,Cc,Dc) in NEGATIVE feedback
    u = K y. Returns the closed-loop A matrix on [x; xc]."""
    Ac, Bc, Cc, Dc = K
    n, nc = A.shape[0], Ac.shape[0]
    Acl = np.block([[A + B @ Dc @ C, B @ Cc],
                    [Bc @ C,          Ac]])
    return Acl


def hinf_norm(A, B, C, D=None, w=None):
    """Peak |G(jw)| on a grid. Grid-based, so it is a lower bound; the grid is
    dense enough here that it tracks the true peak closely."""
    if w is None:
        w = np.geomspace(1e-2, 1e3, 3000)
    n = A.shape[0]
    D = np.zeros((C.shape[0], B.shape[1])) if D is None else D
    return max(np.linalg.norm(C @ np.linalg.solve(1j*wi*np.eye(n) - A, B) + D, 2)
               for wi in w)
