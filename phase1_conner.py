"""
Phase 1, re-anchored to the Conner/Tang/Dowell/Virgin (1997) Duke rig.

Replaces the earlier synthetic 0.30 m / 6 Hz / 2-accelerometer harness. What
the real rig actually does, from the paper:

  SENSORS      three RVDTs, independent: pitch angle at the upper end of the
               pitch axis, plunge from the upper support block, flap rotation
               relative to the main wing from a micro-RVDT on the flap axis.
               DISPLACEMENTS, not accelerations. C = [I3 | 0], D = 0.
               Angular RVDT linearity quoted at 0.26 percent.
  SAMPLING     4000 points per channel at 800 points/s -> 5 s records.
  EXCITATION   B&K 4810 mini-shaker at MIDSPAN, QUARTER CHORD, with a B&K 8200
               force transducer measuring the force at the joint, so the input
               is genuinely measured. Transfer functions 0-50 Hz, 10 averages,
               Hanning. Note a = -0.5 puts the elastic axis AT the quarter
               chord, so the shaker force acts through the EA: it drives plunge
               and produces no direct pitching moment. That is a real property
               of this rig, not a modelling shortcut.
  DAMPING      Conner's measured log-dec ratios, on by default.

Two plant modes:
  linear   -> exact ZOH marching of A(U), fast
  freeplay -> RK4 + Henon switching from freeplay.py

STATE: z = [alpha, beta, h, alphadot, betadot, hdot, x1, x2]
"""
import numpy as np
from scipy.linalg import expm, inv
from scipy.signal import butter, sosfiltfilt, chirp

from conner3dof import ConnerSection, state_space, aero_matrices
from freeplay import integrate as fp_integrate

FS = 800.0            # Conner: 800 points/s
NSAMP = 4000          # Conner: 4000 points per channel -> 5 s


def input_matrix(sec, U):
    """Shaker force at the quarter chord = at the elastic axis, positive down.
    Enters the h equation only; no direct pitching moment."""
    A_non, _, _ = aero_matrices(sec, U, 0.0)
    Minv = inv(sec.Ms + A_non)
    B = np.zeros((8, 1))
    B[3:6, 0] = Minv @ np.array([0.0, 0.0, 1.0])
    return B


def output_matrices():
    """Three RVDTs: alpha [rad], beta [rad], h [m]."""
    C = np.zeros((3, 8)); C[0, 0] = C[1, 1] = C[2, 2] = 1.0
    return C, np.zeros((3, 1))


# ---------------------------------------------------------------------------
def excite_sweep(t, f0=0.2, f1=50.0, amp=1.0):
    """Their GVT excitation: auto-sweeping sinusoid over the 0-50 Hz band."""
    return amp*chirp(t, f0, t[-1], f1, method="linear")


def excite_raps(t, n, rng, width=0.004, amp=1.0):
    u = np.zeros_like(t)
    for tk in np.linspace(t[-1]*0.03, t[-1]*0.97, n):
        u += amp*np.exp(-0.5*((t - tk)/width)**2)
    return u


# ---------------------------------------------------------------------------
def simulate_linear(sec, U, u_func, T=NSAMP/FS, fs=FS, oversample=16,
                    rng=None, snr_db=40.0, adc_bits=16,
                    structural_damping=True):
    """Analog-rate ZOH march, sensor noise at the analog rate, anti-alias,
    decimate to fs. Noise is inflated so snr_db is the ACHIEVED in-band SNR."""
    rng = rng or np.random.default_rng(0)
    A = state_space(sec, U, structural_damping, "stiff")
    B = input_matrix(sec, U)
    C, D = output_matrices()

    fs_hi = fs*oversample
    dt = 1.0/fs_hi
    n = int(round(T*fs_hi))
    t = np.arange(n)*dt
    u = np.atleast_2d(u_func(t))

    M = np.zeros((9, 9)); M[:8, :8] = A; M[:8, 8:] = B
    Md = expm(M*dt); Ad, Bd = Md[:8, :8], Md[:8, 8:]
    x = np.zeros((8, n))
    for i in range(n - 1):
        x[:, i+1] = Ad @ x[:, i] + Bd @ u[:, i]
    y = C @ x + D @ u
    return _sample(t, u, y, fs, oversample, rng, snr_db, adc_bits, U, A, B, C, D)


def simulate_freeplay(sec, U, delta, T=NSAMP/FS, fs=FS, dt=2e-4, rng=None,
                      snr_db=40.0, adc_bits=16, beta0_mult=1.0,
                      settle=12.0, structural_damping=True):
    """LCO record: settle first, then capture a 5 s window at fs."""
    rng = rng or np.random.default_rng(0)
    d = delta
    tt, zz, ns = fp_integrate(sec, U, d, T=settle + T, dt=dt,
                              z0=[0, beta0_mult*d, 0, 0, 0, 0, 0, 0],
                              structural_damping=structural_damping,
                              record_from=settle)
    C, _ = output_matrices()
    y = C @ zz.T
    idx = np.arange(0, len(tt), max(1, int(round(1/(fs*dt)))))
    ts, ys = tt[idx], y[:, idx]
    rms = np.std(ys, axis=1)*10**(-snr_db/20)
    ys = ys + rng.standard_normal(ys.shape)*np.atleast_2d(rms).T
    if adc_bits:
        q = 4*np.max(np.abs(ys))/2**adc_bits
        ys = np.round(ys/q)*q
    return dict(t=ts, u=np.zeros((1, len(ts))), y=ys, U=U, fs=fs,
                delta=d, switches=ns, mode="freeplay")


def _sample(t, u, y, fs, oversample, rng, snr_db, adc_bits, U, A, B, C, D):
    rms = np.std(y, axis=1)*10**(-snr_db/20)
    analog = np.asarray(rms)/np.sqrt(0.8/oversample)
    yn = y + rng.standard_normal(y.shape)*np.atleast_2d(analog).T
    sos = butter(8, 0.8*(fs/2)/(fs*oversample/2), output="sos")
    yf = sosfiltfilt(sos, yn, axis=1)
    idx = np.arange(0, y.shape[1], oversample)
    ys = yf[:, idx]
    if adc_bits:
        q = 4*np.max(np.abs(ys))/2**adc_bits
        ys = np.round(ys/q)*q
    return dict(t=t[idx], u=u[:, idx], y=ys, y_clean=y[:, idx], U=U, fs=fs,
                A=A, B=B, C=C, D=D, mode="linear",
                modes=modal_truth(A))


def modal_truth(A):
    ev = np.linalg.eigvals(A)
    osc = ev[ev.imag > 1e-9]
    osc = osc[np.argsort(osc.imag)]
    return np.array([[abs(s)/2/np.pi, -s.real/abs(s)] for s in osc])


def build_dataset(sec, UF, ratios=(0.30, 0.50, 0.65, 0.80, 0.90, 0.95),
                  mode="sweep", seed=0, snr_db=40.0, n_raps=6):
    rng = np.random.default_rng(seed)
    out = []
    for r in ratios:
        if mode == "sweep":
            f = lambda t: excite_sweep(t)
        else:
            f = lambda t: excite_raps(t, n_raps, rng)
        d = simulate_linear(sec, r*UF, f, rng=rng, snr_db=snr_db)
        d["ratio"] = r
        out.append(d)
    return out
