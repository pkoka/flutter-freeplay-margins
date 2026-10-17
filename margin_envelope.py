"""
Amplitude dependence measured WITHIN a single ring-down, then turned into a
bound rather than a scatter of point estimates.

Part 1 (tightening). A decaying record on a freeplay structure is not a
stationary system: early in the decay the flap slams the stops (effective
stiffness -> K_beta), late in the decay it rattles inside the gap (effective
stiffness -> 0). Sliding a window down one record therefore sweeps the
effective stiffness continuously. This removes the awkward "shaker amplitude in
synthetic units" axis from the earlier experiment and replaces it with a
physical one: flap amplitude relative to the gap.

Part 2 (the bound). The describing-function family K_eff = gamma*K_beta,
gamma in [0,1], contains EVERY quasi-static amplitude the structure can be in.
Minimising the instability speed over that family gives a guaranteed lower
bound on the speed at which the structure can first sustain motion - one
number, not a scatter.
"""
import numpy as np
from scipy.signal import butter, sosfiltfilt
from scipy.optimize import brentq, minimize_scalar
from conner3dof import ConnerSection, state_space, ss_flutter
from phase1_conner import input_matrix, output_matrices
from freeplay import integrate
from sysid import po_moesp, modal_params

FS, DT_SIM = 800.0, 1.25e-4   # dt must divide 1/FS exactly; 2e-4 gave 833.33 Hz


def ringdown(sec, U, delta, amp, T=10.0, snr_db=40.0, seed=0):
    rng = np.random.default_rng(seed)
    B = input_matrix(sec, U); C, _ = output_matrices()
    uf = lambda t: amp*np.exp(-0.5*((t - 0.25)/0.004)**2)
    t, z, ns = integrate(sec, U, delta, T=T, dt=DT_SIM, u_func=uf, B=B)
    step = int(round(1/(FS*DT_SIM)))
    assert abs(1/(step*DT_SIM) - FS) < 1e-9, "dt does not divide 1/FS"
    idx = np.arange(0, len(t), step)
    ts, zz = t[idx], z[idx]
    y = (C @ zz.T)
    y = y + rng.standard_normal(y.shape)*np.atleast_2d(
        np.std(y, axis=1)*10**(-snr_db/20)).T
    sos = butter(6, 50/(FS/2), output="sos")
    return ts, sosfiltfilt(sos, y, axis=1), zz


def windowed_id(ts, y, zz, delta, win=1.5, hop=0.5, order=6, seed=0):
    """Slide a window down the decay; report (rms |beta|*sqrt2/delta, f, zeta).

    The identification is output-only (a free ring-down has no exogenous
    input), so PO-MOESP is fed a deterministic negligible dither rather than an
    exactly zero input, which would make its regression singular. The dither is
    SEEDED: with an unseeded draw the reported modal parameters change from run
    to run, which makes any figure built on this function unreproducible.
    """
    n = int(win*FS); h = int(hop*FS)
    rng = np.random.default_rng(seed)
    dither = 1e-12*rng.standard_normal((1, n))
    out = []
    for s in range(0, len(ts) - n, h):
        sl = slice(s, s + n)
        yb = y[:, sl]
        if np.std(yb) < 1e-12:
            continue
        amp = np.sqrt(np.mean(zz[sl, 1]**2))*np.sqrt(2)/delta
        u = np.zeros((1, n))
        try:
            M = modal_params(po_moesp(u + dither, yb,
                                      i=40, order=order)[0], 1/FS, fmax=40)
        except Exception:
            continue
        if len(M) == 0:
            continue
        k = np.argmin(M[:, 1])
        out.append((amp, M[k, 0], M[k, 1]))
    return np.array(out)


# --- Part 2: the describing-function envelope -----------------------------
def onset_of_gamma(sec, gamma, sd=True, lo=0.2, hi=40.0):
    As = state_space(sec, 1.0, sd, "stiff")   # placeholder, rebuilt per U

    def mx(U):
        A = state_space(sec, U, sd, "free")
        return np.linalg.eigvals(
            A + gamma*(state_space(sec, U, sd, "stiff") - A)).real.max()
    if mx(lo) > 0:
        return lo
    try:
        return brentq(mx, lo, hi)
    except ValueError:
        return hi


def envelope(sec, ngamma=201, sd=True):
    gs = np.linspace(0.0, 1.0, ngamma)
    Us = np.array([onset_of_gamma(sec, g, sd) for g in gs])
    return gs, Us


if __name__ == "__main__":
    sec = ConnerSection()
    UF, _ = ss_flutter(sec, np.linspace(1., 40., 3000), True)
    print(f"true linear U_f = {UF:.3f} m/s\n")

    print("PART 1 - amplitude dependence within ONE ring-down, U = 0.26 U_f")
    for dg in (1.15, 2.12):
        d = np.deg2rad(dg)
        ts, y, zz = ringdown(sec, 0.26*UF, d, 200.0, T=12.0, seed=11)
        W = windowed_id(ts, y, zz, d)
        print(f"  gap {dg} deg:   |beta|/gap -> zeta")
        for a, f, z in W[::2]:
            print(f"      {a:7.3f}   f={f:6.2f} Hz   zeta={z:+.4f}")

    print("\nPART 2 - describing-function envelope over K_eff = gamma*K_beta")
    gs, Us = envelope(sec)
    i = np.argmin(Us)
    print(f"  guaranteed lower bound: U >= {Us[i]:.3f} m/s = {Us[i]/UF:.4f} U_f"
          f"   (at gamma = {gs[i]:.3f})")
    print(f"  upper end (gamma = 1, nominal linear): {Us[-1]/UF:.4f} U_f")
    for g in (0.0, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0):
        print(f"    gamma={g:5.2f}  onset = {onset_of_gamma(sec,g)/UF:.4f} U_f")
