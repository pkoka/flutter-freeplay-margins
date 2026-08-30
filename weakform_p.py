"""
Item 2, done properly: weak-form fitting at an order that can actually
represent the plant.

Item 1 established that a 3-output, 8-state plant needs a THIRD-order matrix
ODE (or d = 3 delays). The memo's weak-form result used a second-order model,
which item 1 proved admits no exact representation. So the derivative-noise
question has to be re-asked at order 3.

MODEL.  y''' + C2 y'' + C1 y' + C0 y = B u        (y in R^3, u in R^1)

WEAK FORM.  Multiply by a compactly supported test function and integrate by
parts i times to move every derivative onto the test function:

    int phi y^(i) dt = (-1)^i int phi^(i) y dt

so the regression needs NO numerical derivative of y at any order, including
the third. The test function is (1-s^2)^q with q > 3 so that phi''' exists and
vanishes at the window edges. Derivatives of phi are taken exactly, as
polynomial coefficients, not numerically.

STRONG FORM.  Same model, same data, y' y'' y''' by finite differences. This is
the comparison that matters: at order 3 the strong form must differentiate
noisy data three times.
"""
import numpy as np
from numpy.polynomial import Polynomial as P
from scipy.linalg import lstsq, eig


def testfun_derivs(nr, q, order):
    """phi(s) = (1-s^2)^q on [-1,1] and its exact derivatives up to `order`."""
    base = (P([1, 0, -1]))**q
    s = np.linspace(-1, 1, 2*nr + 1)
    return [base.deriv(i)(s) if i else base(s) for i in range(order + 1)]


def weak_fit_p(q_meas, u, dt, p=3, nr=24, nwin=600, qpow=9, reg=None):
    """Fit y^(p) + sum_{i<p} C_i y^(i) = B u in weak form. Returns C_0..C_{p-1}, B."""
    nq, N = q_meas.shape
    m = u.shape[0]
    phis = testfun_derivs(nr, qpow, p)
    L = len(phis[0])
    rt = nr*dt
    centers = np.linspace(nr, N - nr - 1, nwin).astype(int)

    ncol = p*nq + m
    A = np.zeros((nq*len(centers), nq*ncol))
    b = np.zeros(nq*len(centers))
    for wi, ci in enumerate(centers):
        sl = slice(ci - nr, ci - nr + L)
        Q, Uw = q_meas[:, sl], u[:, sl]
        # int phi^(i) y dt / rt^i, with the (-1)^i from integration by parts
        I = [((-1)**i)*(Q*phis[i]).sum(1)*dt/rt**i for i in range(p + 1)]
        Iu = (Uw*phis[0]).sum(1)*dt
        feats = np.concatenate([np.concatenate([I[i] for i in range(p)]), Iu])
        for r in range(nq):
            A[wi*nq + r, r*ncol:(r+1)*ncol] = feats
            b[wi*nq + r] = -I[p][r]
    sc = np.sqrt((A**2).sum(0)); sc[sc == 0] = 1.0
    th = (lstsq(A/sc, b, cond=reg)[0]/sc).reshape(nq, ncol)
    Cs = [th[:, i*nq:(i+1)*nq] for i in range(p)]
    return Cs, th[:, p*nq:]


def strong_fit_p(q_meas, u, dt, p=3):
    """Same model, derivatives by finite differences."""
    nq = q_meas.shape[0]
    d = [q_meas]
    for _ in range(p):
        d.append(np.gradient(d[-1], dt, axis=1))
    X = np.vstack([d[i] for i in range(p)] + [u]).T
    th = lstsq(X, -d[p].T)[0].T
    Cs = [th[:, i*nq:(i+1)*nq] for i in range(p)]
    return Cs, th[:, p*nq:]


def modes_p(Cs, fmax=40.0):
    """Modal parameters from the companion form of y^(p) + sum C_i y^(i) = 0."""
    p = len(Cs); n = Cs[0].shape[0]
    A = np.zeros((p*n, p*n))
    A[:-n, n:] = np.eye((p-1)*n)
    A[-n:, :] = -np.hstack(Cs)
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
