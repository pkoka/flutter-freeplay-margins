"""
Path A, item 1: is the freeplay residual sparse in state coordinates?

Item 1 of the original list closed with a hard negative: for the Conner rig an
exact order-3 output-delay model needs 24 of 27 nonzero coefficients, so sparse
regression has nothing to find. That result is about representing the WHOLE
plant in delay coordinates. This script asks the different question the
original proposal actually needed: given the nominal linear model, is the
RESIDUAL sparse in state coordinates?

Parts:
  1a  exact structural identity, r = C_beta * E * sat_delta(beta)
  1b  sparsity count against the delay-coordinate floor
  1c  zero-noise sparse recovery, true delta in the dictionary
  1d  zero-noise sparse recovery, true delta NOT in the dictionary, and
      whether the coefficient-weighted centroid estimates delta

STATE: z = [alpha, beta, h, alphadot, betadot, hdot, x1, x2]
"""
import numpy as np
from scipy.linalg import inv
from conner3dof import ConnerSection, aero_matrices, ss_flutter
from freeplay import regions, integrate

np.set_printoptions(precision=4, suppress=True)


def injection_direction(sec, U):
    """E: how a flap hinge moment enters zdot. Nonzero in the three
    acceleration rows only."""
    A_non, _, _ = aero_matrices(sec, U, 0.0)
    Minv = inv(sec.Ms + A_non)
    E = np.zeros(8)
    E[3:6] = Minv @ np.array([0.0, 1.0, 0.0])
    return E


def record(sec, U, delta, beta0_mult=1.0, T=6.0, settle=1.0, dt=2e-4):
    """Exact states and exact state derivatives from the piecewise fields.
    No finite differences and no measurement noise: this isolates the
    representation question from the estimation question."""
    A_s, A_f, ap, an = regions(sec, U, delta)
    t, z, ns = integrate(sec, U, delta, T=settle + T, dt=dt,
                         z0=[0, beta0_mult*delta, 0, 0, 0, 0, 0, 0],
                         record_from=settle)
    def zd(zz):
        b = zz[1]
        if b > delta:  return A_s @ zz + ap
        if b < -delta: return A_s @ zz + an
        return A_f @ zz
    zdot = np.array([zd(zi) for zi in z])
    R = zdot - z @ A_s.T            # residual against the STIFF nominal model
    return t, z, R, ns, (A_s, A_f, ap, an)


def build_library(zz, dgrid):
    """Linear terms in every state, odd and even polynomials in beta, the
    usual non-smooth candidates, two cross terms, and a dictionary of
    saturations at candidate deadband sizes in degrees."""
    names, cols = [], []
    for i, nm in enumerate(["alpha","beta","h","alphad","betad","hd","x1","x2"]):
        names.append(nm); cols.append(zz[:, i])
    for p in (2, 3, 5, 7):
        names.append(f"beta^{p}"); cols.append(zz[:, 1]**p)
    names += ["|beta|", "sign(beta)", "alpha*beta", "beta*betad"]
    cols += [np.abs(zz[:, 1]), np.sign(zz[:, 1]),
             zz[:, 0]*zz[:, 1], zz[:, 1]*zz[:, 4]]
    for d in dgrid:
        names.append(f"sat{d:.2f}")
        cols.append(np.clip(zz[:, 1], -np.deg2rad(d), np.deg2rad(d)))
    return np.column_stack(cols), names


def stlsq(Th, R, thr=1e-3, rows=(3, 4, 5), n_it=25):
    """STLSQ solved per output row, thresholded against that row's own largest
    coefficient. Rows outside `rows` are left at zero: 1a establishes they are
    identically zero, and thresholding a numerically-zero column against its
    own max keeps pure roundoff."""
    sc = np.linalg.norm(Th, axis=0)
    Tn = Th/sc
    Xi = np.zeros((Th.shape[1], R.shape[1]))
    for k in rows:
        x = np.linalg.lstsq(Tn, R[:, k], rcond=None)[0]
        for _ in range(n_it):
            small = np.abs(x) < thr*np.max(np.abs(x))
            x[small] = 0.0
            big = ~small
            if big.sum():
                x[big] = np.linalg.lstsq(Tn[:, big], R[:, k], rcond=None)[0]
        Xi[:, k] = x
    return Xi/sc[:, None]


def main():
    sec = ConnerSection()
    UF, fF = ss_flutter(sec, np.linspace(1.0, 40.0, 800))
    print(f"stiff-region flutter speed U_f = {UF:.4f} m/s at {fF:.4f} Hz\n")

    DELTA = np.deg2rad(2.12)
    U = 0.73*UF
    E = injection_direction(sec, U)
    t, z, R, ns, (A_s, A_f, ap, an) = record(sec, U, DELTA)

    # ---------------------------------------------------------------- 1a
    print("="*74); print("1a  STRUCTURAL IDENTITY"); print("="*74)
    dA = A_f - A_s
    nz = np.argwhere(np.abs(dA) > 1e-30)
    print(f"  A_free - A_stiff nonzero entries : {len(nz)}  at {[tuple(x) for x in nz]}")
    print(f"  || dA[:,1] - C_beta E ||_inf     : {np.max(np.abs(dA[:,1]-sec.Cb*E)):.3e}")
    print(f"  || a_pos  - C_beta delta E ||_inf: {np.max(np.abs(ap-sec.Cb*DELTA*E)):.3e}")
    print(f"  || a_neg  + C_beta delta E ||_inf: {np.max(np.abs(an+sec.Cb*DELTA*E)):.3e}")
    Rp = np.outer(sec.Cb*np.clip(z[:,1],-DELTA,DELTA), E)
    sv = np.linalg.svd(R, compute_uv=False)
    print(f"  samples {len(t)}, switching events {ns}")
    print(f"  || R - C_beta sat(beta) E ||_inf : {np.max(np.abs(R-Rp)):.3e}   (scale {np.max(np.abs(R)):.3e})")
    print(f"  singular values of R             : {sv[0]:.3e} {sv[1]:.3e} {sv[2]:.3e}")
    print(f"  numerical rank of R              : {np.linalg.matrix_rank(R, tol=1e-9*sv[0])}")
    print(f"  fraction of samples inside the deadband : {100*np.mean(np.abs(z[:,1])<=DELTA):.1f}%\n")

    # ---------------------------------------------------------------- 1b
    print("="*74); print("1b  SPARSITY COUNT"); print("="*74)
    bmax = np.rad2deg(np.abs(z[:,1]).max())
    dg = np.arange(1.0, min(bmax, 3.65), 0.2)
    Th, nm = build_library(z, np.append(dg, 2.12))
    print("  output-delay representation (item 1, closed):")
    print("    27 coefficients, freedom l(dl-n)=3, minimum nonzero 24  ->  89% dense")
    print("  state-coordinate residual representation:")
    print(f"    library {Th.shape[1]} terms, Xi is {Th.shape[1]}x8 = {Th.shape[1]*8} coefficients")
    print(f"    true model: 1 active term, 3 nonzero rows  ->  {100*3/(Th.shape[1]*8):.2f}% dense\n")

    # ---------------------------------------------------------------- 1c
    print("="*74); print("1c  ZERO-NOISE RECOVERY, true delta IN the dictionary"); print("="*74)
    Xi = stlsq(Th, R)
    act = [nm[i] for i in range(len(nm)) if np.any(Xi[i] != 0)]
    i = nm.index("sat2.12")
    rr = np.linalg.norm(R[:,3:6]-(Th@Xi)[:,3:6])/np.linalg.norm(R[:,3:6])
    print(f"  active terms          : {act}")
    print(f"  nonzero coefficients  : {int(np.sum(Xi[:,3:6]!=0))} in rows 3-5")
    print(f"  recovered / (C_b E)   : {Xi[i,3:6]/(sec.Cb*E[3:6])}")
    print(f"  relative residual     : {rr:.3e}\n")

    # ---------------------------------------------------------------- 1d
    print("="*74); print("1d  ZERO-NOISE RECOVERY, true delta NOT in the dictionary"); print("="*74)
    print(f"{'d_true':>7} {'U/Uf':>6} {'b0/d':>6} {'grid':>6} {'lo':>5} {'hi':>5} "
          f"{'centroid':>9} {'err%':>8} {'sum/CbE':>9} {'nnz':>5} {'rel res':>10}")
    for dtrue in (2.12, 1.50, 3.00):
        D = np.deg2rad(dtrue)
        for r, b0, lo, step in ((0.73,1.0,1.00,0.05), (0.73,1.0,0.50,0.05),
                                (0.73,3.0,0.50,0.05), (0.90,1.0,0.50,0.05),
                                (0.60,2.0,0.50,0.05), (0.73,1.0,0.50,0.10)):
            Uu = r*UF
            Ee = injection_direction(sec, Uu)
            _, zz, RR, _, _ = record(sec, Uu, D, beta0_mult=b0)
            bm = np.rad2deg(np.abs(zz[:,1]).max())
            if bm < dtrue*1.05:
                continue
            hi = min(bm, 3.65) if lo == 1.00 else bm
            g = np.arange(lo, hi, step); g = g[np.abs(g-dtrue) > 1e-9]
            TH, NM = build_library(zz, g)
            X = stlsq(TH, RR)
            sats = [(float(NM[k][3:]), X[k,4]) for k in range(len(NM))
                    if NM[k].startswith("sat") and X[k,4] != 0]
            tot = sum(c for _, c in sats)
            cen = sum(d*c for d, c in sats)/tot if tot else np.nan
            res = np.linalg.norm(RR[:,3:6]-(TH@X)[:,3:6])/np.linalg.norm(RR[:,3:6])
            print(f"{dtrue:>7.2f} {r:>6.2f} {b0:>6.1f} {step:>6.2f} {lo:>5.2f} {hi:>5.2f} "
                  f"{cen:>9.3f} {100*(cen-dtrue)/dtrue:>8.2f} {tot/(sec.Cb*Ee[4]):>9.4f} "
                  f"{int(np.sum(X[:,3:6]!=0)):>5d} {res:>10.2e}")
    print()
    print("  VERDICT. The residual is exactly 1-sparse and rank 1, and sparse")
    print("  regression recovers it to machine precision WHEN the true deadband")
    print("  is in the dictionary. When it is not, the fit spreads across the")
    print("  whole saturation dictionary, the coefficient sum is not conserved,")
    print("  and the coefficient-weighted centroid moves with the dictionary")
    print("  range and the operating point. Delta must be measured, not fitted.")


if __name__ == "__main__":
    main()
