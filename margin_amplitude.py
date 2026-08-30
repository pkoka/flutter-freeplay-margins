"""
Does the identified flutter margin depend on the amplitude you tested at?

Standard subcritical flutter testing: excite at several speeds below the
boundary, identify modal damping, extrapolate damping to zero, call that the
flutter speed. On a linear structure the answer is a property of the structure.
On a structure with control-surface freeplay it cannot be, because the
effective control-surface stiffness the test "sees" depends on how far the flap
swings relative to the gap. Small raps stay inside the gap and identify the
FREE system; large raps spend most of the cycle against the stops and identify
the STIFF system. Those two have different flutter speeds.

This script measures that spread on the Conner rig, using the same shaker,
sensors, sampling and identification as everything else in Phase 2.
"""
import numpy as np
from scipy.signal import butter, sosfiltfilt
from conner3dof import ConnerSection, state_space, ss_flutter
from phase1_conner import input_matrix, output_matrices
from freeplay import integrate
from sysid import po_moesp, modal_params

DT_SIM, FS = 2e-4, 800.0


def rap_record(sec, U, delta, amp, T=6.0, settle=0.0, snr_db=40.0, seed=0):
    """Single shaker rap on the freeplay plant; return the sampled ring-down."""
    rng = np.random.default_rng(seed)
    B = input_matrix(sec, U)
    C, _ = output_matrices()
    uf = lambda t: amp*np.exp(-0.5*((t - 0.25)/0.004)**2)
    t, z, ns = integrate(sec, U, delta, T=T, dt=DT_SIM, u_func=uf, B=B,
                         record_from=settle)
    step = int(round(1/(FS*DT_SIM)))
    idx = np.arange(0, len(t), step)
    ts, y = t[idx], (C @ z.T)[:, idx]
    u = np.array([[uf(tt) for tt in ts]])
    rms = np.std(y, axis=1)*10**(-snr_db/20)
    y = y + rng.standard_normal(y.shape)*np.atleast_2d(rms).T
    sos = butter(6, 50/(FS/2), output="sos")
    return u, sosfiltfilt(sos, y, axis=1), ns, np.abs(z[:, 1]).max()/delta


def identify_damping(u, y, order=6):
    """Damping of the mode that goes unstable (the flutter mode)."""
    try:
        M = modal_params(po_moesp(u, y, i=40, order=order)[0], 1/FS, fmax=40)
    except Exception:
        return np.nan, np.nan
    if len(M) == 0:
        return np.nan, np.nan
    k = np.argmin(M[:, 1])          # least-damped mode
    return M[k, 0], M[k, 1]


def extrapolate_Uf(Us, zetas):
    """Linear extrapolation of damping to zero: the standard subcritical
    flutter-test practice. Uses the last few points, as flight test does."""
    m = np.isfinite(zetas)
    Us, zetas = np.asarray(Us)[m], np.asarray(zetas)[m]
    if len(Us) < 3:
        return np.nan
    p = np.polyfit(Us[-4:], zetas[-4:], 1)
    if p[0] >= 0:
        return np.inf                # damping still rising: no crossing implied
    return -p[1]/p[0]


if __name__ == "__main__":
    sec = ConnerSection()
    UF, _ = ss_flutter(sec, np.linspace(1., 40., 3000), True)
    delta = np.deg2rad(2.12)
    speeds = np.array([0.10, 0.14, 0.18, 0.22, 0.26])*UF     # below LCO onset
    # amplitudes chosen to span the gap: the shaker acts through the
    # elastic axis, so it takes ~10 N/m of impulse to move the flap one
    # gap width. Below that the flap never leaves region (2).
    amps = [1.0, 5.0, 10.0, 20.0, 50.0, 200.0, 1000.0]

    print(f"Conner rig, gap = 2.12 deg, true linear U_f = {UF:.2f} m/s")
    print(f"free-region instability at 0.286 U_f; tests run at "
          f"{speeds[0]/UF:.2f}-{speeds[-1]/UF:.2f} U_f, all below LCO onset\n")
    print("  rap amp   max|beta|/delta   f [Hz]   zeta at highest speed   "
          "extrapolated U_f      as U/U_f")
    for amp in amps:
        zs, fs_, bmax = [], [], []
        for U in speeds:
            u, y, ns, bm = rap_record(sec, U, delta, amp, seed=11)
            f, z = identify_damping(u, y)
            zs.append(z); fs_.append(f); bmax.append(bm)
        Uf_hat = extrapolate_Uf(speeds, zs)
        r = Uf_hat/UF if np.isfinite(Uf_hat) else np.inf
        print(f"  {amp:8.1e}   {np.max(bmax):9.2f}      {fs_[-1]:6.2f}   "
              f"{zs[-1]:+18.4f}   {Uf_hat:14.2f}   {r:9.3f}")
