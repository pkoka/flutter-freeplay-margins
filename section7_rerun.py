"""Re-run of memo Section 7 with the corrected mode matcher and at an order
that can represent the plant. Same 12 cases as the original: three airspeeds
crossed with four noise levels."""
import numpy as np
from scipy.linalg import expm
from scipy.signal import butter, sosfiltfilt
from conner3dof import ConnerSection, ss_flutter, state_space
from phase1_conner import simulate_linear, excite_sweep, input_matrix, output_matrices
from weakform import weak_fit, weak_modes
from weakform_p import weak_fit_p, strong_fit_p, modes_p
from sysid import po_moesp, modal_params, era_from_decay
from scoring import score, truth_modes

FS, DT = 800.0, 1/800.
RATIOS = (0.65, 0.90, 0.95)
SNRS = (60., 40., 30., 20.)
SOS = butter(6, 50/(FS/2), output="sos")
# tuned on a held-out seed, per SNR
W2 = {60.: (10, 5), 40.: (10, 9), 30.: (13, 9), 20.: (32, 3)}
W3 = {60.: (40, 5), 40.: (40, 5), 30.: (32, 5), 20.: (24, 5)}


def strong_fit_p2(q, u, dt):
    qd = np.gradient(q, dt, axis=1); qdd = np.gradient(qd, dt, axis=1)
    from scipy.linalg import lstsq
    X = np.vstack([qd, q, u]).T
    th = lstsq(X, -qdd.T)[0].T
    return th[:, :3], th[:, 3:6]


def decay_record(sec, U, snr, seed=5):
    A = state_space(sec, U, True, "stiff"); B = input_matrix(sec, U)
    C, _ = output_matrices(); n = 4000; Ad = expm(A*DT)
    x = np.zeros((8, n)); x[:, 0] = B[:, 0]/np.linalg.norm(B)
    for k in range(n-1):
        x[:, k+1] = Ad @ x[:, k]
    y = C @ x
    rng = np.random.default_rng(seed)
    y = y + rng.standard_normal(y.shape)*np.atleast_2d(np.std(y, axis=1)*10**(-snr/20)).T
    return sosfiltfilt(SOS, y, axis=1)


def run():
    sec = ConnerSection()
    UF, _ = ss_flutter(sec, np.linspace(1., 40., 3000), True)
    res = {k: {s: [] for s in SNRS} for k in
           ("moesp8", "moesp6", "era", "w2", "s2", "w3", "s3")}
    for snr in SNRS:
        for r in RATIOS:
            T = truth_modes(sec, r*UF)
            d0 = simulate_linear(sec, r*UF, lambda t: excite_sweep(t),
                                 rng=np.random.default_rng(3), snr_db=snr)
            y = sosfiltfilt(SOS, d0["y"], axis=1)
            u = sosfiltfilt(SOS, d0["u"], axis=1)
            f = lambda M: score(M, T)[1]
            res["moesp8"][snr].append(f(modal_params(po_moesp(u, y, i=40, order=8)[0], DT, fmax=40)))
            res["moesp6"][snr].append(f(modal_params(po_moesp(u, y, i=40, order=6)[0], DT, fmax=40)))
            res["era"][snr].append(f(era_from_decay(decay_record(sec, r*UF, snr), DT,
                                                    order=8, nmarkov=1200, r=400, s=780)))
            nr, qp = W2[snr]
            _, C2, K2, _ = weak_fit(y, u, DT, nr=nr, nwin=600, p=qp)
            res["w2"][snr].append(f(weak_modes(C2, K2)))
            C2s, K2s = strong_fit_p2(y, u, DT)
            res["s2"][snr].append(f(weak_modes(C2s, K2s)))
            nr, qp = W3[snr]
            Cs, _ = weak_fit_p(y, u, DT, p=3, nr=nr, nwin=800, qpow=qp)
            res["w3"][snr].append(f(modes_p(Cs)))
            Cs, _ = strong_fit_p(y, u, DT, p=3)
            res["s3"][snr].append(f(modes_p(Cs)))
    return res, UF


if __name__ == "__main__":
    res, UF = run()
    names = [("moesp8", "PO-MOESP order 8"), ("moesp6", "PO-MOESP order 6"),
             ("era", "ERA, free decay"), ("w3", "weak form, order 3"),
             ("s3", "strong form, order 3"), ("w2", "weak form, order 2"),
             ("s2", "strong form, order 2")]
    print(f"U_f = {UF:.3f} m/s.  Mean damping error %, "
          f"averaged over U/U_f = {RATIOS}\n")
    print(f"  {'method':22s}" + "".join(f"{s:>9.0f} dB" for s in SNRS) + "     all")
    for k, lab in names:
        row = [np.nanmean(res[k][s]) for s in SNRS]
        allm = np.nanmean([v for s in SNRS for v in res[k][s]])
        cells = "".join(f"{v:12.2f}" if np.isfinite(v) else f"{'n/a':>12s}" for v in row)
        print(f"  {lab:22s}{cells}{allm:8.2f}")
