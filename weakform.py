"""
Weak-form least squares identification of a second-order model.

This is the WSINDy idea applied honestly to this problem. The plant is
    M qddot + C qdot + K q = B u ,   q = [alpha, beta, h]
and we measure only q (three RVDTs). Multiplying by a compactly supported test
function phi and integrating by parts moves BOTH derivatives onto phi:

    int phi*M*qddot dt = int phi'' *M*q dt
    int phi*C*qdot  dt = -int phi' *C*q dt

so no numerical differentiation of noisy data is required anywhere. That is
the actual selling point of the weak form, and it is why it belongs in this
comparison rather than a nonlinear library: the plant IS linear inside each
freeplay region, so a nonlinear library would only add variance.

The model is deliberately misspecified in one known way: a second-order model
in [alpha,beta,h] cannot represent the two Wagner aerodynamic lag states. That
unmodelled content is exactly what should show up as identification error, and
in Phase 3 as uncertainty.
"""
import numpy as np
from scipy.linalg import lstsq, eig


def _testfun(nr, p=9):
    """phi(s) = (1-s^2)^p on s in [-1,1], with its first two derivatives."""
    s = np.linspace(-1, 1, 2*nr + 1)
    phi = (1 - s**2)**p
    d1 = p*(1 - s**2)**(p-1)*(-2*s)
    d2 = p*(p-1)*(1 - s**2)**(p-2)*(4*s**2) + p*(1 - s**2)**(p-1)*(-2)
    return phi, d1, d2


def weak_fit(q, u, dt, nr=24, nwin=400, p=9, reg=None):
    """NOTE on nr: the test function acts as a band-pass whose usable range is
    roughly f <~ sqrt(p)/(2 pi nr dt). With nr=120 at 800 Hz that is ~1 Hz,
    which low-passes away modes at 5-19 Hz entirely. Size nr to the HIGHEST
    mode you need, not to the record length."""
    """q:(3,N) displacements, u:(m,N) input. Returns (Mn, Cn, Kn, Bn) with M
    normalized to identity (only the ratios are identifiable from q and u)."""
    nq, N = q.shape
    m = u.shape[0]
    phi, d1, d2 = _testfun(nr, p)
    L = len(phi)
    rt = nr*dt                              # window half-width in seconds
    centers = np.linspace(nr, N - nr - 1, nwin).astype(int)

    # each window gives nq equations:  qdd_w + C q_w' + K q_w = B u_w
    # in weak form with M = I after normalization:
    #   int phi'' q  +  (-int phi' q) C^T + (int phi q) K^T = (int phi u) B^T
    A = np.zeros((nq*len(centers), nq*(2*nq + m)))
    b = np.zeros(nq*len(centers))
    for wi, cidx in enumerate(centers):
        sl = slice(cidx - nr, cidx - nr + L)
        Q = q[:, sl]; Uw = u[:, sl]
        # phi is defined on s = (t-tc)/rt, so d/dt = (d/ds)/rt. The chain-rule
        # factor is the WINDOW HALF-WIDTH rt = nr*dt, not dt.
        I2 = (Q*d2).sum(1)*dt/rt**2             # int phi'' q
        I1 = -(Q*d1).sum(1)*dt/rt               # -int phi' q
        I0 = (Q*phi).sum(1)*dt
        Iu = (Uw*phi).sum(1)*dt
        feats = np.concatenate([I1, I0, Iu])    # (2nq+m,)
        for r in range(nq):
            A[wi*nq + r, r*(2*nq+m):(r+1)*(2*nq+m)] = feats
            b[wi*nq + r] = -I2[r]
    # column scaling: the C-block, K-block and B-block features differ by
    # ~3 orders of magnitude, giving cond ~ 1e5 unscaled. Solve the balanced
    # problem and undo the scaling afterwards.
    sc = np.sqrt((A**2).sum(0)); sc[sc == 0] = 1.0
    th = (lstsq(A/sc, b, cond=reg)[0]/sc).reshape(nq, 2*nq + m)
    Cn, Kn, Bn = th[:, :nq], th[:, nq:2*nq], th[:, 2*nq:]
    return np.eye(nq), Cn, Kn, Bn


def weak_modes(Cn, Kn, fmax=60.0):
    """Modal frequencies [Hz] and damping from the normalized (I, C, K)."""
    n = Cn.shape[0]
    A = np.block([[np.zeros((n, n)), np.eye(n)], [-Kn, -Cn]])
    lam = eig(A)[0]
    out = []
    for s in lam:
        if s.imag <= 1e-9:
            continue
        f = abs(s)/(2*np.pi)
        if f > fmax:
            continue
        out.append((f, -s.real/abs(s)))
    out.sort()
    return np.array(out) if out else np.zeros((0, 2))
