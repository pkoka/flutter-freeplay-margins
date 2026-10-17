"""
Stage 2 endpoint. Does a data-measured uncertainty weight change the
controller you get, and does it change whether that controller survives the
true plant?

PLANT. Free-region Conner model, which is the one that goes unstable, at
0.286 U_f nominal. Control input is a commanded flap hinge moment, the same
channel the freeplay enters. Measurements are the three RVDTs.

DESIGN. Mixed sensitivity with unit control weight, synthesised on the
NOMINAL model. Robust stability to input multiplicative uncertainty is then
evaluated as  ||W T_i||_inf < 1  with T_i the complementary sensitivity at the
plant input.

TEST. Run both controllers on the TRUE perturbed plant WITH freeplay, which
switches between A_free and A_stiff and was never in any design model.
"""
import numpy as np
from conner3dof import ConnerSection, state_space, ss_flutter, aero_matrices
from phase1_conner import output_matrices
from scipy.linalg import inv, expm
from hinf import hinf_syn, lft_closed_loop
from stage2a_harness import perturbed_section

W_RAD = np.geomspace(1e-2, 1e3, 1500)


KA = 100.0    # actuator authority. B scaled to unit norm alone is far too weak
              # against this plant and inflates gamma from ~2.6 to ~179 at U_f.


def plant(sec, U, region="free", ka=None):
    A = state_space(sec, U, True, region)
    A_non, _, _ = aero_matrices(sec, U, 0.0)
    B = np.zeros((8, 1)); B[3:6, 0] = inv(sec.Ms + A_non) @ np.array([0., 1., 0.])
    B = B/np.linalg.norm(B)*(KA if ka is None else ka)
    C, _ = output_matrices()
    C = np.diag([1/0.02, 1/0.05, 1/0.002]) @ C   # normalise rad, rad, m
    return A, B, C


def genplant(A, B, C, wp):
    B1 = np.hstack([B, np.zeros((8, 3))])
    C1 = np.vstack([wp*C, np.zeros((1, 8))])
    D12 = np.vstack([np.zeros((3, 1)), np.ones((1, 1))])
    D21 = np.hstack([np.zeros((3, 1)), np.eye(3)])
    return A, B1, B, C1, C, D12, D21


def Ti_norm(A, B, C, K, Wv):
    """||W(w) T_i(jw)||_inf, T_i = complementary sensitivity at the PLANT INPUT."""
    Ac, Bc, Cc, _ = K
    n, nc = A.shape[0], Ac.shape[0]
    Acl = lft_closed_loop(A, B, C, K)
    if np.linalg.eigvals(Acl).real.max() >= 0:
        return np.inf
    Bcl = np.vstack([B, np.zeros((nc, 1))])          # disturbance at plant input
    Ccl = np.hstack([np.zeros((1, n)), Cc])          # controller output
    peak = 0.0
    for wi, Wi in zip(W_RAD, Wv):
        T = (Ccl @ np.linalg.solve(1j*wi*np.eye(n+nc) - Acl, Bcl)).item()
        peak = max(peak, Wi*abs(T))
    return peak


def design(sec, U, wp):
    A, B, C = plant(sec, U)
    K, g = hinf_syn(*genplant(A, B, C, wp))   # (A,B1,B2,C1,C2,D12,D21)
    return A, B, C, K, g


if __name__ == "__main__":
    sec_n = ConnerSection(); sec_t, _ = perturbed_section()
    UF, _ = ss_flutter(sec_n, np.linspace(1., 40., 800))
    Uopen = 0.286*UF
    print(f"U_f (stiff, nominal) {UF:.3f} m/s;  free-region open-loop instability "
          f"at {Uopen:.3f} m/s = 0.286 U_f\n")

    print("H-infinity synthesis on the nominal free-region plant")
    print(f"{'U/Uf':>6} {'wp':>6} {'gamma':>9} {'cl stable':>10} {'||T_i||':>9}")
    ones = np.ones_like(W_RAD)
    for r in (0.30, 0.40, 0.60, 0.80, 1.00):
        for wp in (0.3, 1.0, 3.0):
            A, B, C, K, g = design(sec_n, r*UF, wp)
            if K is None:
                print(f"{r:>6.2f} {wp:>6.1f} {'-':>9} {'no design':>10} {'-':>9}")
                continue
            st = np.linalg.eigvals(lft_closed_loop(A, B, C, K)).real.max()
            print(f"{r:>6.2f} {wp:>6.1f} {g:>9.3f} "
                  f"{'yes' if st < 0 else 'NO':>10} {Ti_norm(A,B,C,K,ones):>9.3f}")


# ---------------------------------------------------------------------------
def closed_loop_freeplay(sec, U, delta, K, Kscale, T=8.0, dt=1.25e-4, z0=None):
    """Controller in the loop with the true switching plant. Fixed-step RK4,
    no Henon: this test is about stability, not about LCO amplitude to four
    digits, and the switching error at this step size is far below the
    growth/decay rates being measured."""
    Ac, Bc, Cc, _ = K
    nc = Ac.shape[0]
    A_s = state_space(sec, U, True, "stiff")
    A_f = state_space(sec, U, True, "free")
    A_non, _, _ = aero_matrices(sec, U, 0.0)
    Minv = inv(sec.Ms + A_non)
    off = Minv @ np.array([0.0, sec.Cb*delta, 0.0])
    ap = np.zeros(8); ap[3:6] = off
    an = -ap
    B = np.zeros((8, 1)); B[3:6, 0] = Minv @ np.array([0., 1., 0.])
    B = B/np.linalg.norm(B)*KA
    C, _ = output_matrices()
    C = np.diag([1/0.02, 1/0.05, 1/0.002]) @ C
    x = np.zeros(8 + nc)
    x[:8] = z0 if z0 is not None else np.array([0, 2*delta, 0, 0, 0, 0, 0, 0])

    def f(x):
        z, xc = x[:8], x[8:]
        b = z[1]
        A, a = (A_s, ap) if b > delta else (A_s, an) if b < -delta else (A_f, np.zeros(8))
        y = C @ z
        u = Cc @ xc
        return np.concatenate([A @ z + a + (B @ u), Ac @ xc + Bc @ y])

    n = int(T/dt)
    peak = np.zeros(n//100 + 1)
    for i in range(n):
        k1 = f(x); k2 = f(x + dt/2*k1); k3 = f(x + dt/2*k2); k4 = f(x + dt*k3)
        x = x + dt/6*(k1 + 2*k2 + 2*k3 + k4)
        if not np.all(np.isfinite(x)) or np.abs(x[:8]).max() > 1e4:
            return np.inf, i*dt
        if i % 100 == 0:
            peak[i//100] = np.abs(x[1])/delta
    tail = peak[-len(peak)//4:]
    return float(tail.max()), T


# ---------------------------------------------------------------------------
def attractor_trace(sec, U, delta, K, T=8.0, dt=1e-4, b0=2.0, nrec=400):
    """Closed-loop amplitude envelope over time. A limit cycle is an attractor:
    constant amplitude, independent of the initial condition. Slow decay is
    neither. Run this from at least two initial conditions before calling any
    closed-loop result stable, because the closed loop here is BISTABLE.

    Stepping is by matrix exponential, exact within each region. Explicit RK4
    cannot be used: near-optimal H-infinity central controllers carry poles
    around 1e5 rad/s here, so RK4 needs dt < 3e-5 and blows up at any usable
    step size. The exponential step is unconditionally stable and one matvec.
    """
    Ac, Bc, Cc, _ = K
    nc = Ac.shape[0]
    n = 8 + nc
    A_non, _, _ = aero_matrices(sec, U, 0.0)
    Minv = inv(sec.Ms + A_non)
    B = np.zeros((8, 1)); B[3:6, 0] = Minv @ np.array([0., 1., 0.])
    B = B/np.linalg.norm(B)*KA
    C = np.diag([1/0.02, 1/0.05, 1/0.002]) @ output_matrices()[0]
    off = np.zeros(8); off[3:6] = Minv @ np.array([0., sec.Cb*delta, 0.])

    def build(region):
        A = state_space(sec, U, True, region)
        return np.block([[A + np.zeros((8, 8)), B @ Cc], [Bc @ C, Ac]])

    steps = {}
    for tag, region, aff in (("free", "free", np.zeros(8)),
                             ("pos", "stiff", off), ("neg", "stiff", -off)):
        Acl = build(region)
        Phi = expm(Acl*dt)
        b = np.concatenate([aff, np.zeros(nc)])
        # exact affine part, solved rather than inverted (Acl can be near singular)
        try:
            v = np.linalg.solve(Acl, (Phi - np.eye(n)) @ b)
        except np.linalg.LinAlgError:
            v = dt*b
        steps[tag] = (Phi, v)

    x = np.zeros(n); x[1] = b0*delta
    nstep = int(T/dt); env = []; blk = max(1, nstep//nrec); cur = 0.0
    for i in range(nstep):
        b = x[1]
        Phi, v = steps["pos"] if b > delta else steps["neg"] if b < -delta else steps["free"]
        x = Phi @ x + v
        if not np.all(np.isfinite(x)) or np.abs(x[:8]).max() > 1e4:
            return np.array(env), False
        cur = max(cur, abs(x[1])/delta)
        if i % blk == blk - 1:
            env.append(cur); cur = 0.0
    return np.array(env), True
