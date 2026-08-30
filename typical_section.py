"""
Phase 0: 2-DOF pitch-plunge typical section with Theodorsen unsteady aero.

Three independent flutter solvers on the same plant:
  1. k-method (V-g), exact Theodorsen C(k)
  2. p-k method, exact Theodorsen C(k)
  3. time-domain state space with R.T. Jones 2-lag Wagner approximation, eig(A(U))

Solvers 1 and 2 are frequency domain and assume harmonic motion.
Solver 3 makes no harmonic assumption. Agreement between all three is the
Phase 0 pass criterion.

CONVENTIONS (fix these in your head, do not mix them anywhere downstream):
  h      plunge, positive DOWN, units of length
  theta  pitch, positive NOSE UP, radians
  b      semichord
  a      elastic axis offset from MIDCHORD, in semichords, positive AFT
  e      mass center offset from MIDCHORD, in semichords, positive AFT
  x_th   static unbalance = e - a, in semichords
  r_th   radius of gyration about the EA, in semichords
  mu     mass ratio = m / (pi rho b^2)
  sigma  omega_h / omega_theta
  k      reduced frequency = omega b / U
  s_nd   nondimensional time = U t / b, in semichords traveled

Theodorsen lift and moment about the EA:
  L   = pi rho b^2 (hdd + U thd - b a thdd)
        + 2 pi rho U b C(k) [hd + U th + b(1/2 - a) thd]
  M_ea= pi rho b^2 [b a hdd - U b(1/2 - a) thd - b^2(1/8 + a^2) thdd]
        + 2 pi rho U b^2 (a + 1/2) C(k) [hd + U th + b(1/2 - a) thd]

Structural EOM (h positive down so lift enters the plunge row with a + sign):
  m hdd + m b x_th thdd + k_h h + L = 0
  m b x_th hdd + I_th thdd + k_th th - M_ea = 0
"""

import numpy as np
from scipy.special import hankel2
from scipy.linalg import eig, inv

# R.T. Jones approximation of the Wagner function
# phi(s) = 1 - psi1 exp(-eps1 s) - psi2 exp(-eps2 s)
PSI1, EPS1 = 0.165, 0.0455
PSI2, EPS2 = 0.335, 0.3000


# ----------------------------------------------------------------------------
# plant definition
# ----------------------------------------------------------------------------
class TypicalSection:
    """2-DOF typical section. Everything downstream consumes this object."""

    def __init__(self, a=-0.2, e=-0.1, mu=20.0, r2=6/25, sigma=0.4,
                 b=1.0, rho=1.0, omega_th=1.0):
        self.a, self.e, self.mu, self.r2, self.sigma = a, e, mu, r2, sigma
        self.b, self.rho, self.omega_th = b, rho, omega_th

        self.x_th = e - a                      # static unbalance, semichords
        self.m = mu * np.pi * rho * b**2       # mass per unit span
        self.I_th = self.m * b**2 * r2         # inertia about the EA
        self.omega_h = sigma * omega_th
        self.k_h = self.m * self.omega_h**2
        self.k_th = self.I_th * omega_th**2

        # structural matrices, q = [h, theta]
        self.Ms = np.array([[self.m, self.m * b * self.x_th],
                            [self.m * b * self.x_th, self.I_th]])
        self.Ks = np.diag([self.k_h, self.k_th])

        # aerodynamic added mass, same for every solver
        self.Ma = np.pi * rho * np.array([[b**2, -b**3 * a],
                                          [-b**3 * a, b**4 * (1/8 + a**2)]])
        self.Mtot = self.Ms + self.Ma

    def __repr__(self):
        return (f"TypicalSection(a={self.a}, e={self.e}, mu={self.mu}, "
                f"r2={self.r2}, sigma={self.sigma})")


def theodorsen(k):
    """Exact Theodorsen C(k) from Hankel functions of the second kind."""
    k = np.asarray(k, dtype=float)
    k = np.where(k < 1e-10, 1e-10, k)          # C(0) = 1, avoid the pole
    H1, H0 = hankel2(1, k), hankel2(0, k)
    return H1 / (H1 + 1j * H0)


def theodorsen_jones(k):
    """C(k) implied by the Jones 2-lag Wagner fit. Use to isolate fit error."""
    p = 1j * np.asarray(k, dtype=float)
    return 1.0 - PSI1 * p / (p + EPS1) - PSI2 * p / (p + EPS2)


# ----------------------------------------------------------------------------
# frequency-domain aero matrices
# ----------------------------------------------------------------------------
def aero_matrices(ts, C):
    """Aero damping and stiffness with the airspeed factored out.

    Returns Chat, Khat such that the full aero contribution to the LHS is
    Ma * s^2 + U * Chat * s + U^2 * Khat, with C frozen at C(k).
    """
    b, rho, a = ts.b, ts.rho, ts.a
    pr = np.pi * rho

    Chat = np.array([
        [2*pr*b*C,
         pr*b**2 + 2*pr*b**2*(0.5 - a)*C],
        [-2*pr*b**2*(a + 0.5)*C,
         pr*b**3*(0.5 - a) - 2*pr*b**3*(a + 0.5)*(0.5 - a)*C]
    ], dtype=complex)

    Khat = np.array([
        [0.0, 2*pr*b*C],
        [0.0, -2*pr*b**2*(a + 0.5)*C]
    ], dtype=complex)

    return Chat, Khat


# ----------------------------------------------------------------------------
# solver 1: k-method (V-g)
# ----------------------------------------------------------------------------
def k_method(ts, k_vals, C_func=theodorsen):
    """Classical V-g. Returns arrays (nk, 2) of U, omega, g per branch."""
    b = ts.b
    nk = len(k_vals)
    U = np.full((nk, 2), np.nan)
    om = np.full((nk, 2), np.nan)
    g = np.full((nk, 2), np.nan)

    for i, k in enumerate(k_vals):
        C = complex(C_func(k))
        Chat, Khat = aero_matrices(ts, C)
        # Ks (1 + i g) q = omega^2 [Mtot - i (b/k) Chat - (b/k)^2 Khat] q
        B = ts.Mtot - 1j*(b/k)*Chat - (b/k)**2 * Khat
        lam, _ = eig(inv(ts.Ks) @ B)           # lam = (1 + i g) / omega^2

        good = lam.real > 1e-12
        w = np.where(good, 1.0/np.sqrt(np.abs(lam.real)), np.nan)
        gg = np.where(good, lam.imag/lam.real, np.nan)

        order = np.argsort(w)                  # branch 0 = lower frequency
        om[i] = w[order]
        g[i] = gg[order]
        U[i] = b*om[i]/k

    return U, om, g


# ----------------------------------------------------------------------------
# solver 2: p-k method
# ----------------------------------------------------------------------------
def _quad_eig(M, Cm, K):
    """Eigenvalues of M s^2 + Cm s + K = 0 via the companion linearization."""
    n = M.shape[0]
    A = np.block([[np.zeros((n, n)), np.eye(n)],
                  [-inv(M) @ K, -inv(M) @ Cm]])
    return eig(A)[0]


def pk_method(ts, U_vals, C_func=theodorsen, tol=1e-10, itmax=200):
    """p-k. Returns (nU, 2) arrays of omega and damping ratio zeta."""
    b = ts.b
    nU = len(U_vals)
    om = np.full((nU, 2), np.nan)
    zeta = np.full((nU, 2), np.nan)

    # seed with the in-vacuo frequencies
    k_seed = [ts.omega_h*b/U_vals[0], ts.omega_th*b/U_vals[0]]

    for i, U in enumerate(U_vals):
        for mode in (0, 1):
            k = max(k_seed[mode]*U_vals[max(i-1, 0)]/U, 1e-6) if i else k_seed[mode]
            for _ in range(itmax):
                C = complex(C_func(k))
                Chat, Khat = aero_matrices(ts, C)
                roots = _quad_eig(ts.Mtot, U*Chat, ts.Ks + U**2 * Khat)
                pos = roots[roots.imag > 1e-9]
                if len(pos) < 2:
                    k = np.nan
                    break
                pos = pos[np.argsort(pos.imag)]
                s = pos[mode]
                k_new = s.imag*b/U
                if abs(k_new - k) < tol*max(abs(k), 1.0):
                    k = k_new
                    break
                k = 0.5*k + 0.5*k_new          # damped update, keeps it stable
            if np.isnan(k):
                continue
            C = complex(C_func(k))
            Chat, Khat = aero_matrices(ts, C)
            roots = _quad_eig(ts.Mtot, U*Chat, ts.Ks + U**2 * Khat)
            pos = roots[roots.imag > 1e-9]
            pos = pos[np.argsort(pos.imag)]
            s = pos[mode]
            om[i, mode] = abs(s)
            zeta[i, mode] = -s.real/abs(s)
            k_seed[mode] = k

    return om, zeta


# ----------------------------------------------------------------------------
# solver 3: time-domain state space, Jones lags
# ----------------------------------------------------------------------------
def state_space(ts, U):
    """A(U) for z = [h, theta, hdot, thetadot, x1, x2]. This is the object
    every later phase consumes, so keep the state ordering fixed."""
    b, rho, a = ts.b, ts.rho, ts.a
    pr = np.pi * rho

    # structural + aero stiffness acting on [h, theta]
    S = np.array([
        [ts.k_h,  pr*U**2*b],
        [0.0,     ts.k_th - pr*U**2*b**2*(a + 0.5)]
    ])
    # damping acting on [hdot, thetadot]
    D = np.array([
        [pr*U*b,                 pr*U*b**2*(1.5 - a)],
        [-pr*U*b**2*(a + 0.5),   pr*U*b**3*(0.5 - a)**2]
    ])
    # lag feed-through acting on [x1, x2]
    G = np.array([
        [2*pr*U**2*PSI1*EPS1,               2*pr*U**2*PSI2*EPS2],
        [-2*pr*U**2*b*(a + 0.5)*PSI1*EPS1,  -2*pr*U**2*b*(a + 0.5)*PSI2*EPS2]
    ])

    Minv = inv(ts.Mtot)
    A = np.zeros((6, 6))
    A[0:2, 2:4] = np.eye(2)
    A[2:4, 0:2] = -Minv @ S
    A[2:4, 2:4] = -Minv @ D
    A[2:4, 4:6] = -Minv @ G
    # xdot_i = -eps_i (U/b) x_i + hdot + U theta + b(1/2 - a) thetadot
    for j, eps in enumerate((EPS1, EPS2)):
        A[4 + j, 1] = U
        A[4 + j, 2] = 1.0
        A[4 + j, 3] = b*(0.5 - a)
        A[4 + j, 4 + j] = -eps*U/b
    return A


def ss_root_locus(ts, U_vals):
    """Eigenvalues of A(U) over the sweep. Returns (nU, 6) complex array."""
    return np.array([np.sort_complex(eig(state_space(ts, U))[0]) for U in U_vals])


# ----------------------------------------------------------------------------
# flutter point extraction
# ----------------------------------------------------------------------------
def _cross(x, y):
    """First upward zero crossing of y(x), linearly interpolated."""
    for i in range(len(x) - 1):
        if np.isnan(y[i]) or np.isnan(y[i+1]):
            continue
        if y[i] < 0.0 <= y[i+1]:
            f = -y[i]/(y[i+1] - y[i])
            return x[i] + f*(x[i+1] - x[i]), f, i
    return np.nan, np.nan, None


def flutter_k(ts, k_vals, C_func=theodorsen):
    """Flutter point from the V-g curves. Sweeps k downward so U increases."""
    U, om, g = k_method(ts, k_vals, C_func)
    best = (np.inf, np.nan)
    for j in (0, 1):
        Uc, f, i = _cross(U[:, j], g[:, j])
        if not np.isnan(Uc) and Uc < best[0]:
            best = (Uc, om[i, j] + f*(om[i+1, j] - om[i, j]))
    return best


def flutter_pk(ts, U_vals, C_func=theodorsen):
    om, zeta = pk_method(ts, U_vals, C_func)
    best = (np.inf, np.nan)
    for j in (0, 1):
        Uc, f, i = _cross(U_vals, -zeta[:, j])
        if not np.isnan(Uc) and Uc < best[0]:
            best = (Uc, om[i, j] + f*(om[i+1, j] - om[i, j]))
    return best


def flutter_ss(ts, U_vals):
    """Flutter = first U where any eigenvalue of A(U) crosses into the RHP."""
    mx = np.array([eig(state_space(ts, U))[0].real.max() for U in U_vals])
    Uc, f, i = _cross(U_vals, mx)
    if np.isnan(Uc):
        return np.inf, np.nan
    ev = eig(state_space(ts, Uc))[0]
    wf = abs(ev[np.argmax(ev.real)])
    return Uc, wf


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    # Hodges & Pierce flutter Example 1 (Georgia Tech aeroelasticity resources):
    # a = -1/5, e = -1/10, mu = 20, r = 2/5, sigma = 2/5
    ts = TypicalSection(a=-0.2, e=-0.1, mu=20.0, r2=6/25, sigma=0.4)
    print(ts)
    print(f"  x_theta = {ts.x_th:.4f}   omega_h/omega_theta = {ts.sigma}")
    print(f"  in-vacuo coupled freqs = "
          f"{np.sort(np.sqrt(np.linalg.eigvals(inv(ts.Ms) @ ts.Ks).real))}")

    k_vals = np.geomspace(3.0, 0.01, 4000)
    U_vals = np.linspace(0.05, 6.0, 3000)

    Uk, wk = flutter_k(ts, k_vals, theodorsen)
    Up, wp = flutter_pk(ts, U_vals, theodorsen)
    Upj, wpj = flutter_pk(ts, U_vals, theodorsen_jones)
    Us, ws = flutter_ss(ts, U_vals)

    b, wt = ts.b, ts.omega_th
    print("\n  method                              V* = U/(b w_th)   w_f/w_th")
    print(f"  k-method,  exact C(k)                    {Uk/(b*wt):8.4f}   {wk/wt:8.4f}")
    print(f"  p-k,       exact C(k)                    {Up/(b*wt):8.4f}   {wp/wt:8.4f}")
    print(f"  p-k,       Jones C(k)                    {Upj/(b*wt):8.4f}   {wpj/wt:8.4f}")
    print(f"  state space, Jones lags, eig(A(U))       {Us/(b*wt):8.4f}   {ws/wt:8.4f}")

    print(f"\n  p-k exact vs p-k Jones   : {100*abs(Upj-Up)/Up:6.3f} %  (Jones fit error)")
    print(f"  p-k Jones vs state space : {100*abs(Us-Upj)/Upj:6.3f} %  (harmonic assumption)")
    print(f"  k-method vs p-k exact    : {100*abs(Uk-Up)/Up:6.3f} %")
