"""
Mode matching for identification scoring.

The greedy nearest-frequency matcher used earlier is not safe. Identified
models routinely contain spurious heavily damped poles, and greedy matching
will assign one of those to a true mode and report a meaningless error. That
produced a 787% score for PO-MOESP at 60 dB in a run where it had in fact
recovered every mode to four digits.

Two changes, applied identically to every method:
  1. Drop candidates with zeta > ZMAX. The true modes of this plant have
     zeta <= 0.21, so a pole at zeta = 0.99 is a numerical artifact, not a
     structural mode. This is a physical filter, not a per-method tuning knob.
  2. Optimal assignment (Hungarian) on relative frequency distance, instead of
     greedy nearest. Greedy can consume a candidate that another true mode
     needed.
"""
import numpy as np
from scipy.optimize import linear_sum_assignment
from conner3dof import state_space

ZMAX = 0.5


def match_modes(M, T, zmax=ZMAX):
    """Return the rows of M assigned to each row of T, or None if impossible."""
    if len(M) == 0:
        return None
    M = M[M[:, 1] < zmax]
    if len(M) < len(T):
        return None
    cost = np.abs(M[:, 0][None, :] - T[:, 0][:, None])/T[:, 0][:, None]
    r, c = linear_sum_assignment(cost)
    return M[c]


def score(M, T, zmax=ZMAX):
    """Mean relative error in frequency and damping, in percent."""
    Mm = match_modes(M, T, zmax)
    if Mm is None:
        return np.nan, np.nan
    return (np.mean(np.abs(Mm[:, 0]/T[:, 0] - 1))*100,
            np.mean(np.abs(Mm[:, 1]/T[:, 1] - 1))*100)


def truth_modes(sec, U, fmax=40.0, region="stiff", sd=True):
    """Ground-truth modal frequencies [Hz] and damping from eig A(U)."""
    ev = np.linalg.eigvals(state_space(sec, U, sd, region))
    osc = ev[ev.imag > 1e-9]
    osc = osc[np.argsort(osc.imag)]
    out = [(abs(s)/2/np.pi, -s.real/abs(s)) for s in osc
           if abs(s)/2/np.pi <= fmax]
    return np.array(out)
