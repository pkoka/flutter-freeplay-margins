"""
Rigorous lower bounds on the speed at which a freeplay structure can sustain
motion, via absolute stability theory.

SETUP. The only difference between the stiff and free regions is the flap
restoring moment, so the plant can be written exactly as a linear system in
feedback with a memoryless nonlinearity:

    zdot = A_free z + B_phi * (-phi(beta)),     beta = C_beta z
    A_stiff = A_free - B_phi * K_beta * C_beta        (verified numerically)

and the freeplay restoring moment obeys, for EVERY beta,

    0 <= phi(beta)/beta <= K_beta                      (sector [0, K_beta])

because the deadband is flat in the middle and has slope K_beta outside it.
This is exact, not an approximation - no averaging, no assumed waveform.

BOUNDS.
  Circle criterion, sector [0,k], G Hurwitz:  Re[G(jw)] > -1/k for all w.
  Popov criterion (valid because freeplay is memoryless and time-invariant,
  which is strictly stronger information than the circle criterion uses):
      exists eta >= 0 :  Re[(1 + j w eta) G(jw)] + 1/k > 0 for all w.
Either condition certifies global asymptotic stability for ALL nonlinearities
in the sector. Below the certified speed no disturbance of any size can produce
sustained motion. That is what a describing function cannot claim.
"""
import numpy as np
from scipy.linalg import inv
from scipy.optimize import brentq, minimize_scalar
from conner3dof import ConnerSection, state_space, aero_matrices


def loop(sec, U, sd=True):
    """Return A_free, B_phi, C_beta, k with A_stiff = A_free - B k C verified."""
    A_free = state_space(sec, U, sd, "free")
    A_stiff = state_space(sec, U, sd, "stiff")
    A_non, _, _ = aero_matrices(sec, U, 0.0)
    Minv = inv(sec.Ms + A_non)
    B = np.zeros((8, 1)); B[3:6, 0] = Minv @ np.array([0.0, 1.0, 0.0])
    C = np.zeros((1, 8)); C[0, 1] = 1.0
    k = sec.Cb
    err = np.abs(A_stiff - (A_free - B*k @ C)).max()
    return A_free, B, C, k, err


def G(A, B, C, w):
    n = A.shape[0]
    return np.array([complex((C @ np.linalg.solve(1j*wi*np.eye(n) - A, B)).item())
                     for wi in w])


def circle_margin(sec, U, w, sd=True):
    """min_w Re[G] + 1/k. Positive => circle criterion certifies stability.
    Returns -inf if A_free is not Hurwitz (criterion needs it)."""
    A, B, C, k, _ = loop(sec, U, sd)
    if np.linalg.eigvals(A).real.max() >= 0:
        return -np.inf
    return float(np.min(G(A, B, C, w).real) + 1/k)


def popov_margin(sec, U, w, sd=True):
    """max over eta >= 0 of min_w Re[(1+j w eta)G] + 1/k."""
    A, B, C, k, _ = loop(sec, U, sd)
    if np.linalg.eigvals(A).real.max() >= 0:
        return -np.inf
    g = G(A, B, C, w)
    def m(eta):
        return float(np.min(g.real - eta*w*g.imag) + 1/k)
    # eta must be searched over MANY decades. Where Im[G] <= 0 for all w the
    # margin grows without bound in eta (certification is easy); where Im[G]
    # goes positive, large eta destroys it. A narrow search range silently
    # returns the range endpoint instead of the optimum - that was a bug.
    etas = np.concatenate([[0.0], np.geomspace(1e-6, 1e8, 300)])
    vals = np.array([m(e) for e in etas])
    best = float(vals.max())
    if vals.argmax() not in (0, len(etas)-1):
        j = vals.argmax()
        r = minimize_scalar(lambda e: -m(e),
                            bounds=(etas[j-1], etas[j+1]), method="bounded")
        best = max(best, -float(r.fun))
    return best


def certified_speed(sec, margin_fn, w, lo=0.2, hi=None, sd=True):
    hi = hi or 6.8
    if margin_fn(sec, hi, w, sd) > 0:
        return hi
    if margin_fn(sec, lo, w, sd) <= 0:
        return np.nan
    return brentq(lambda U: margin_fn(sec, U, w, sd), lo, hi, xtol=1e-4)


if __name__ == "__main__":
    sec = ConnerSection()
    UF = 23.970
    w = np.concatenate([np.geomspace(1e-3, 1e3, 4000)])
    _, _, _, k, err = loop(sec, 10.0)
    print(f"feedback decomposition check: max|A_stiff - (A_free - B k C)| = {err:.2e}")
    print(f"sector upper bound k = K_beta = {k:.5f} N m/rad per unit span\n")

    Uc = certified_speed(sec, circle_margin, w)
    Up = certified_speed(sec, popov_margin, w)
    print(f"  circle criterion certifies stability up to  U = {Uc:.3f} m/s "
          f"= {Uc/UF:.4f} U_f")
    print(f"  Popov  criterion certifies stability up to  U = {Up:.3f} m/s "
          f"= {Up/UF:.4f} U_f")
    print()
    print("  describing-function estimate (heuristic)      0.1612 U_f")
    print("  earliest observed LCO, gap 2.12 deg           0.18   U_f")
    print("  earliest observed LCO, gap 1.15 deg           0.20   U_f")
    print("  Conner numerical / experimental onset         0.18 / 0.23 U_f")
