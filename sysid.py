"""
Phase 2 baselines: ERA (with OKID) and PO-MOESP subspace identification.

These are what the aeroelastic flight-test and modal-analysis communities
actually use. Anything new has to be compared against them on the same data,
scored against the same stored ground truth, or the comparison is meaningless.

Both return a discrete-time (A,B,C,D); modal_params() converts to continuous
frequencies and damping ratios via the matrix logarithm.
"""
import numpy as np
from scipy.linalg import svd, lstsq, qr, logm, eig


# ---------------------------------------------------------------------------
def modal_params(Ad, dt, fmax=None):
    """Continuous-time modal frequencies [Hz] and damping ratios from Ad."""
    lam = eig(Ad)[0]
    s = np.log(lam.astype(complex))/dt
    out = []
    for si in s:
        if si.imag <= 1e-9:
            continue
        f = abs(si)/(2*np.pi)
        z = -si.real/abs(si)
        if fmax is not None and f > fmax:
            continue
        out.append((f, z))
    out.sort()
    return np.array(out) if out else np.zeros((0, 2))


def _hankel_blocks(x, r, s, start=0):
    """Block Hankel of the (l x m) sequence x[k], r block rows, s block cols."""
    l, m = x[0].shape
    H = np.zeros((r*l, s*m))
    for i in range(r):
        for j in range(s):
            H[i*l:(i+1)*l, j*m:(j+1)*m] = x[start + i + j]
    return H


# ---------------------------------------------------------------------------
# ERA
# ---------------------------------------------------------------------------
def era(markov, order, r=None, s=None):
    """Eigensystem Realization Algorithm on Markov parameters h[0..], where
    h[0] = D and h[k] = C A^(k-1) B."""
    D = markov[0]
    h = markov[1:]
    n_h = len(h)
    r = r or n_h//2
    s = s or n_h - r
    H0 = _hankel_blocks(h, r, s, 0)
    H1 = _hankel_blocks(h, r, s, 1)
    U, S, Vt = svd(H0, full_matrices=False)
    n = order
    Un, Sn, Vn = U[:, :n], S[:n], Vt[:n, :].T
    Sh = np.diag(1/np.sqrt(Sn))
    Ad = Sh @ Un.T @ H1 @ Vn @ Sh
    l, m = h[0].shape
    Obs = Un @ np.diag(np.sqrt(Sn))
    Ctr = np.diag(np.sqrt(Sn)) @ Vn.T
    Cd, Bd = Obs[:l, :], Ctr[:, :m]
    return Ad, Bd, Cd, D, S


def okid(u, y, nm=None, reg=1e-8):
    """Observer/Kalman filter ID: recover system Markov parameters from a
    general (non-impulsive) input. u: (m,N), y: (l,N)."""
    m, N = u.shape
    l = y.shape[0]
    nm = nm or min(200, N//5)
    # regression matrix V of [u_k ; (u,y)_{k-1..k-nm}]
    rows = m + nm*(m + l)
    V = np.zeros((rows, N))
    V[:m] = u
    for i in range(1, nm + 1):
        blk = m + (i-1)*(m+l)
        V[blk:blk+m, i:] = u[:, :N-i]
        V[blk+m:blk+m+l, i:] = y[:, :N-i]
    G = lstsq(V.T, y.T, cond=reg)[0].T          # observer Markov params
    D = G[:, :m]
    Yb = [G[:, m + i*(m+l): m + i*(m+l) + m] for i in range(nm)]      # Ybar1
    Yb2 = [G[:, m + i*(m+l) + m: m + (i+1)*(m+l)] for i in range(nm)]  # Ybar2
    # recursion h_k = Ybar1_k + sum_{i=1}^{k} Ybar2_i h_{k-i}
    hs = [D]
    for k in range(nm):
        hk = Yb[k].copy()
        for i in range(k + 1):
            hk = hk + Yb2[i] @ hs[k - i]
        hs.append(hk)
    return hs


# ---------------------------------------------------------------------------
# PO-MOESP subspace identification
# ---------------------------------------------------------------------------
def po_moesp(u, y, i, order):
    """Past-output MOESP. u:(m,N), y:(l,N). i = block rows."""
    m, N = u.shape
    l = y.shape[0]
    j = N - 2*i + 1
    def blk(x, i0, nb):
        d = x.shape[0]
        H = np.zeros((nb*d, j))
        for k in range(nb):
            H[k*d:(k+1)*d, :] = x[:, i0 + k: i0 + k + j]
        return H
    Up, Uf = blk(u, 0, i), blk(u, i, i)
    Yp, Yf = blk(y, 0, i), blk(y, i, i)
    Zp = np.vstack([Up, Yp])
    # LQ of [Uf; Zp; Yf]
    L = qr(np.vstack([Uf, Zp, Yf]).T, mode="economic")[1].T
    a, b = i*m, i*m + Zp.shape[0]
    L32 = L[b:, a:b]
    U_, S_, _ = svd(L32, full_matrices=False)
    n = order
    Obs = U_[:, :n] @ np.diag(np.sqrt(S_[:n]))
    Cd = Obs[:l, :]
    Ad = lstsq(Obs[:-l, :], Obs[l:, :])[0]
    # B, D by linear least squares on the simulation error
    def sim(B, D):
        x = np.zeros(n); out = np.zeros((l, N))
        for k in range(N):
            out[:, k] = Cd @ x + D @ u[:, k]
            x = Ad @ x + B @ u[:, k]
        return out
    nb = n*m + l*m
    Phi = np.zeros((l*N, nb))
    base = sim(np.zeros((n, m)), np.zeros((l, m)))
    for q in range(nb):
        e = np.zeros(nb); e[q] = 1.0
        B = e[:n*m].reshape(n, m); D = e[n*m:].reshape(l, m)
        Phi[:, q] = (sim(B, D) - base).ravel()
    th = lstsq(Phi, y.ravel() - base.ravel())[0]
    return Ad, th[:n*m].reshape(n, m), Cd, th[n*m:].reshape(l, m), S_


# ---------------------------------------------------------------------------
def era_from_decay(y, dt, order=8, nmarkov=600, r=200, s=380, fmax=60.0):
    """ERA run directly on a free-decay (ring-down) record, which is how
    flutter flight test actually uses it: rap the structure, watch it ring
    down, realize from the response. Pairing ERA with OKID on a swept-sine
    record instead is what made it look terrible - that was my error, not
    ERA's. Band-limit y before calling."""
    h = [np.zeros((y.shape[0], 1))] + [y[:, k:k+1] for k in range(nmarkov)]
    Ad = era(h, order=order, r=r, s=s)[0]
    return modal_params(Ad, dt, fmax=fmax)
