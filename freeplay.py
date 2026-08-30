"""
Piecewise-linear freeplay integration for the Conner 3-DOF section.

Conner Fig. 2: the control-surface restoring moment is
    region (1)  beta < -delta :  M = K_beta (beta + delta)
    region (2)  |beta| < delta:  M = 0
    region (3)  beta > +delta :  M = K_beta (beta - delta)
so regions (1) and (3) share the SAME A matrix (full K_beta) and differ only
in the constant offset vector a; region (2) has K_beta = 0 and a = 0. That is
exactly what Conner states, and the state_space() switch was verified to touch
only the beta column of the three acceleration rows.

    zdot = A_r z + a_r

Switching points are located with Henon's method (Henon 1982), as Conner does:
when a time step crosses beta = +/- delta, swap the roles of the dependent
variable beta and the independent variable t and integrate "in space" from the
current beta back to the exact switching value in a single step. Time is
carried along as an extra integrated variable so the leftover time in the step
is known exactly. Bisection or interpolation on adaptive steps is the inferior
alternative Conner explicitly rejects.

STATE: z = [alpha, beta, h, alphadot, betadot, hdot, x1, x2]   (beta = z[1])
"""
import numpy as np
from scipy.linalg import inv
from conner3dof import ConnerSection, state_space, aero_matrices

IB = 1        # index of beta
IBD = 4       # index of betadot


def regions(sec, U, delta, structural_damping=True):
    """Return (A_stiff, A_free, a_pos, a_neg): the two coefficient matrices
    and the offset vectors for regions (3) and (1)."""
    A_stiff = state_space(sec, U, structural_damping, "stiff")
    A_free = state_space(sec, U, structural_damping, "free")
    # offset: M_beta = C_beta*beta -/+ C_beta*delta. A_stiff already carries
    # C_beta*beta, so the constant part enters the acceleration rows as
    # +/- Minv @ [0, C_beta*delta, 0].
    A_non, _, _ = aero_matrices(sec, U, 0.0)
    Minv = inv(sec.Ms + A_non)
    off = Minv @ np.array([0.0, sec.Cb*delta, 0.0])
    a_pos = np.zeros(8); a_pos[3:6] = off       # region (3), beta > +delta
    a_neg = np.zeros(8); a_neg[3:6] = -off      # region (1), beta < -delta
    return A_stiff, A_free, a_pos, a_neg


def _region_of(beta, delta):
    if beta > delta:
        return 3
    if beta < -delta:
        return 1
    return 2


def _fields(r, A_stiff, A_free, a_pos, a_neg):
    if r == 3:
        return A_stiff, a_pos
    if r == 1:
        return A_stiff, a_neg
    return A_free, np.zeros(8)


def _rk4_t(z, A, a, dt):
    f = lambda y: A @ y + a
    k1 = f(z); k2 = f(z + dt/2*k1); k3 = f(z + dt/2*k2); k4 = f(z + dt*k3)
    return z + dt/6*(k1 + 2*k2 + 2*k3 + k4)


def _rk4_beta(z, t, A, a, dbeta):
    """Henon step: integrate d[z,t]/dbeta = [Az+a, 1]/betadot by dbeta."""
    def f(y):
        zz, bd = y[:8], y[IBD]
        if abs(bd) < 1e-14:
            return None
        d = np.empty(9)
        d[:8] = (A @ zz + a)/bd
        d[8] = 1.0/bd
        return d
    y = np.concatenate([z, [t]])
    k1 = f(y)
    if k1 is None:
        return None
    k2 = f(y + dbeta/2*k1)
    k3 = f(y + dbeta/2*k2)
    k4 = f(y + dbeta*k3)
    if k2 is None or k3 is None or k4 is None:
        return None
    y = y + dbeta/6*(k1 + 2*k2 + 2*k3 + k4)
    return y


def integrate(sec, U, delta, T, dt=2e-4, z0=None, structural_damping=True,
              record_from=0.0, u_func=None, B=None):
    """Time-march with exact switching. Returns t, z (recorded portion).

    u_func(t) -> scalar shaker force, with input matrix B (8,1). The forcing is
    held constant across each RK4 step and across the Henon sub-step, which is
    consistent to the order of the scheme for a smooth input."""
    A_stiff, A_free, a_pos, a_neg = regions(sec, U, delta, structural_damping)
    z = np.zeros(8) if z0 is None else np.array(z0, float)
    t = 0.0
    ts, zs = [], []
    nswitch = 0
    while t < T:
        r = _region_of(z[IB], delta)
        A, a = _fields(r, A_stiff, A_free, a_pos, a_neg)
        if u_func is not None:
            a = a + B[:, 0]*u_func(t)
        zt = _rk4_t(z, A, a, dt)
        rt = _region_of(zt[IB], delta)
        if rt != r:
            # crossed a switching surface: Henon back onto it exactly
            # pick the surface we actually crossed, NOT the sign of beta
            # after the step: leaving region (3) downward still targets
            # +delta, and getting this from zt[IB] aims at the wrong stop.
            target = delta if (r == 3 or rt == 3) else -delta
            y = _rk4_beta(z, t, A, a, target - z[IB])
            if y is not None:
                z_s, t_s = y[:8], y[8]
                z_s[IB] = target                       # enforce the surface
                rem = dt - (t_s - t)
                if 0.0 < rem <= dt:
                    r2 = 3 if (target > 0 and z_s[IBD] > 0) else \
                         1 if (target < 0 and z_s[IBD] < 0) else 2
                    A2, a2 = _fields(r2, A_stiff, A_free, a_pos, a_neg)
                    if u_func is not None:
                        a2 = a2 + B[:, 0]*u_func(t_s)
                    z, t = _rk4_t(z_s, A2, a2, rem), t_s + rem
                else:
                    z, t = z_s, t_s
                nswitch += 1
            else:
                z, t = zt, t + dt
        else:
            z, t = zt, t + dt
        if t >= record_from:
            ts.append(t); zs.append(z.copy())
    return np.array(ts), np.array(zs), nswitch


def rms_amplitudes(t, z, sec, delta):
    """Conner's normalization: pitch and flap by the freeplay size, plunge by
    delta*b. Mean removed, since the LCO can sit off-centre."""
    al, be, h = z[:, 0], z[:, 1], z[:, 2]
    A = np.sqrt(np.mean((al - al.mean())**2))/delta
    B = np.sqrt(np.mean((be - be.mean())**2))/delta
    H = np.sqrt(np.mean((h - h.mean())**2))/(delta*sec.b)
    return A, B, H


def dominant_freq(t, y):
    y = y - y.mean()
    if np.std(y) < 1e-12:
        return np.nan
    n = len(y)
    f = np.fft.rfftfreq(n, t[1] - t[0])
    P = np.abs(np.fft.rfft(y*np.hanning(n)))
    return float(f[1:][np.argmax(P[1:])])
