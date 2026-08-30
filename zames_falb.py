"""
Zames-Falb multipliers: a tighter rigorous bound by using more of what we
actually know about freeplay.

Popov only used "sector [0,k], memoryless, time-invariant". But the freeplay
restoring moment is also:
    ODD                 M(-beta) = -M(beta)   (symmetric gap, no preload)
    MONOTONE            dM/dbeta >= 0 everywhere
    SLOPE-RESTRICTED    dM/dbeta in [0, K_beta]
Zames-Falb multipliers are exactly the class that exploits monotonicity, so
they should - in principle - certify a higher speed than Popov.

CRITERION. For phi monotone, odd, in sector [0,k], the loop is stable if there
exists a multiplier M(s) = 1 - Hhat(s), with h in L1 and ||h||_1 <= 1, s.t.

    Re[ M(jw) (1/k + G(jw)) ] > 0   for all w.

For an ODD nonlinearity the sign restriction on h is dropped, so both causal
and anticausal basis terms may enter with either sign.

PARAMETRIZATION. Basis  psi_i(s) = (rho/(s+rho))^i   (causal)
                        psi_i(s) = (rho/(-s+rho))^i  (anticausal)
each of which has impulse response with unit L1 norm, so ||h||_1 <= sum|c_i|
and the constraint sum|c_i| <= 1 is sufficient. The certification condition is
AFFINE in c, so feasibility is a linear program on a frequency grid.
"""
import numpy as np
from scipy.optimize import linprog, brentq
from conner3dof import ConnerSection
from absolute_stability import loop, G


def basis(w, rho, N):
    """Causal and anticausal first-order-lag powers evaluated on jw."""
    s = 1j*w
    out = []
    for i in range(1, N + 1):
        out.append((rho/(s + rho))**i)
        out.append((rho/(-s + rho))**i)
    return np.array(out)                       # (2N, nw)


def zf_margin(sec, U, w, rhos=(6., 20., 60., 200.), N=4, sd=True):
    """Max over multipliers of the worst-case certification margin.
    > 0 means Zames-Falb certifies absolute stability at this speed."""
    A, B, C, k, _ = loop(sec, U, sd)
    if np.linalg.eigvals(A).real.max() >= 0:
        return -np.inf
    F = 1/k + G(A, B, C, w)
    best = -np.inf
    for rho in rhos:
        P = basis(w, rho, N)                   # (m, nw)
        m = P.shape[0]
        # maximize t  s.t.  Re[F] - sum_i c_i Re[psi_i F] >= t ,  sum|c_i| <= 1
        # variables: [t, c+ (m), c- (m)],  c = c+ - c-
        R = (P*F).real                         # (m, nw)
        nw = len(w)
        Aub = np.zeros((nw + 1, 1 + 2*m))
        Aub[:nw, 0] = 1.0
        Aub[:nw, 1:1+m] = R.T
        Aub[:nw, 1+m:] = -R.T
        bub = np.concatenate([F.real, [1.0]])
        Aub[nw, 1:] = 1.0
        cobj = np.zeros(1 + 2*m); cobj[0] = -1.0
        bnds = [(None, None)] + [(0, None)]*(2*m)
        r = linprog(cobj, A_ub=Aub, b_ub=bub, bounds=bnds, method="highs")
        if r.status == 0:
            best = max(best, float(r.x[0]))
    return best


if __name__ == "__main__":
    sec = ConnerSection(); UF = 23.970
    w = np.geomspace(1e-3, 1e3, 1200)
    print("Zames-Falb margin vs speed (positive = certified):")
    for r in (0.10, 0.14, 0.16, 0.18, 0.20, 0.24, 0.27):
        print(f"   U/U_f={r:.2f}   margin = {zf_margin(sec, r*UF, w):+.5f}")
    lo, hi = 0.05*UF, 0.284*UF
    if zf_margin(sec, hi, w) > 0:
        Uz = hi
    else:
        Uz = brentq(lambda U: zf_margin(sec, U, w), lo, hi, xtol=1e-3)
    print(f"\n  ZAMES-FALB certified speed: {Uz:.3f} m/s = {Uz/UF:.4f} U_f")
    print(f"  Popov (previous rigorous bound):            0.16115 U_f")
    print(f"  earliest observed LCO (2.12 deg gap):       0.18    U_f")


# ---------------------------------------------------------------------------
def combined_margin(sec, U, w, rhos=(6., 20., 60., 200.), N=4,
                    etas=None, sd=True):
    """Combined Popov x Zames-Falb multiplier,  M(s) = (1 + eta s)(1 - Hhat(s)).

    Zames-Falb ALONE fails on this plant, and the reason is structural: with
    ||h||_1 <= 1 the multiplier obeys |M(jw)| <= 2, so it can only supply
    bounded phase advance. Popov's (1 + j w eta) is unbounded, and that is what
    this loop needs. The combined class contains both as special cases (eta = 0
    gives ZF, c = 0 gives Popov) so it can never do worse than either.

    For fixed eta the condition is still affine in c, so it stays an LP; eta is
    gridded outside. Search eta over many decades - a narrow range silently
    returns the endpoint, which is the bug that bit the plain Popov search.
    """
    A, B, C, k, _ = loop(sec, U, sd)
    if np.linalg.eigvals(A).real.max() >= 0:
        return -np.inf
    F = 1/k + G(A, B, C, w)
    etas = etas if etas is not None else np.concatenate(
        [[0.0], np.geomspace(1e-5, 1e6, 45)])
    nw = len(w)
    best = -np.inf
    for eta in etas:
        Fe = (1 + 1j*w*eta)*F
        for rho in rhos:
            P = basis(w, rho, N)
            m = P.shape[0]
            R = (P*Fe).real
            Aub = np.zeros((nw + 1, 1 + 2*m))
            Aub[:nw, 0] = 1.0
            Aub[:nw, 1:1+m] = R.T
            Aub[:nw, 1+m:] = -R.T
            bub = np.concatenate([Fe.real, [1.0]])
            Aub[nw, 1:] = 1.0
            cobj = np.zeros(1 + 2*m); cobj[0] = -1.0
            r = linprog(cobj, A_ub=Aub, b_ub=bub,
                        bounds=[(None, None)] + [(0, None)]*(2*m),
                        method="highs")
            if r.status == 0:
                best = max(best, float(r.x[0]))
    return best
