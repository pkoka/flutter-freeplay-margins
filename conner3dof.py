"""
3-DOF Theodorsen typical section with aileron, built on the Conner/Tang/
Dowell/Virgin (1997) Duke wind tunnel rig.

Aerodynamics: Theodorsen NACA 496 equations (A),(B),(C), taken in the recast
matrix form of NASA/TM-2017-219667 (Perry), eqns (3b)-(5d).
  A_struct, C_struct  structural mass and stiffness
  A_non               apparent mass
  B_non, B_circ       aero damping, noncirculatory and circulatory
  C_non, C_circ       aero stiffness, noncirculatory and circulatory
Theodorsen's normalization divides the alpha and beta equations by M b^2 and
the h equation by M b, so multiplying rows by diag(Mb^2, Mb^2, Mb) recovers
dimensional equations. Verified: that scaling turns A_struct into
[[I_a, I_b+b(c-a)S_b, S_a],[.., I_b, S_b],[S_a, S_b, M]] exactly.

q = [alpha, beta, h],  h positive DOWN, alpha nose up, beta TE down.
"""
import numpy as np
from scipy.special import hankel2
from scipy.linalg import eig, inv


def theodorsen(k):
    k = max(float(k), 1e-10)
    H1, H0 = hankel2(1, k), hankel2(0, k)
    return H1/(H1 + 1j*H0)


def T_constants(a, c):
    ac, s = np.arccos(c), np.sqrt(1 - c**2)
    T = {}
    T[1] = -(2 + c**2)/3*s + c*ac
    T[3] = -(1 - c**2)/8*(5*c**2 + 4) + c*(7 + 2*c**2)*s*ac/4 - (1/8 + c**2)*ac**2
    T[4] = c*s - ac
    T[5] = -(1 - c**2) - ac**2 + 2*c*s*ac
    T[7] = c*(7 + 2*c**2)*s/8 - (1/8 + c**2)*ac
    T[10] = s + ac
    T[11] = (2 - c)*s + (1 - 2*c)*ac
    T[12] = (2 + c)*s - (1 + 2*c)*ac
    T['p'] = -(1 - c**2)**1.5/3        # NACA 496 App. I: p = -(1/3)(1-c^2)^(3/2)
    return T


class ConnerSection:
    """Table A.1 of Conner et al. 1997, all per unit span."""

    def __init__(self, plunge_support_mass=True):
        self.b, self.span = 0.127, 0.52
        self.a, self.c = -0.5, 0.5
        self.M = 1.558                            # wing-aileron mass/length
        self.Sa, self.Sb = 0.08587, 0.00395
        self.Ia, self.Ib = 0.01347, 0.0003264
        self.Ka, self.Kb, self.Kh = 1486., 155., 1809.
        self.za, self.zb, self.zh = 0.01626, 0.0115, 0.0113
        self.rho = 1.225                          # from kappa = pi rho b^2/M
        b, M = self.b, self.M

        # support blocks translate in plunge only -> h-h entry only
        self.Mh = M + (2*0.47485/self.span if plunge_support_mass else 0.0)

        self.Ca, self.Cb, self.Ch = self.Ka*M*b**2, self.Kb*M*b**2, self.Kh*M
        self.Ms = np.array([
            [self.Ia,                        self.Ib + b*(self.c-self.a)*self.Sb, self.Sa],
            [self.Ib + b*(self.c-self.a)*self.Sb, self.Ib,                        self.Sb],
            [self.Sa,                        self.Sb,                             self.Mh]])
        self.Ks = np.diag([self.Ca, self.Cb, self.Ch])
        self.kappa = np.pi*self.rho*b**2/M
        self.T = T_constants(self.a, self.c)

    def modal(self):
        w2, Phi = eig(self.Ks, self.Ms)
        idx = np.argsort(w2.real)
        w = np.sqrt(w2.real[idx]); Phi = Phi.real[:, idx]
        return w, Phi

    def Cs_damping(self):
        """Conner Sec. 3: measured damping ratios treated as modal values,
        assembled as diag(2 m_i w_i zeta_i) and transformed back to [a,b,h]."""
        w, Phi = self.modal()
        mmod = np.diag(Phi.T @ self.Ms @ Phi)
        # zeta ordering: modes come out plunge, pitch, flap by frequency
        zeta = np.array([self.zh, self.za, self.zb])
        Bmod = np.diag(2*mmod*w*zeta)
        Pi = inv(Phi)
        return Pi.T @ Bmod @ Pi


def aero_matrices(sec, U, C):
    """Dimensional aero A, B, K from TM-2017-219667 eqns (3c),(4c),(4d),(5c),(5d).
    Returns matrices already row-scaled by diag(Mb^2, Mb^2, Mb)."""
    a, c, b, kap, T = sec.a, sec.c, sec.b, sec.kappa, sec.T
    pi = np.pi
    scale = np.diag([sec.M*b**2, sec.M*b**2, sec.M*b])

    A_non = kap*np.array([
        [1/8 + a**2,                   -(c-a)*T[1]/pi - T[7]/pi,  -a/b],
        [-(c-a)*T[1]/pi - T[7]/pi,     -T[3]/pi**2,               -T[1]/(pi*b)],
        [-a,                           -T[1]/pi,                   1/b]])

    B_non = kap*(U/b)*np.array([
        [0.5 - a,                          (1/pi)*(-2*T['p'] - (0.5-a)*T[4]), 0.],
        [(1/pi)*(T['p'] - T[1] - T[4]/2),  -T[4]*T[11]/(2*pi**2),             0.],
        [1.,                               -T[4]/pi,                          0.]])

    B_cir = C*kap*(U/b)*np.array([
        [2*(a**2 - 0.25),          -(1/pi)*(a+0.5)*T[11],   -2*(a+0.5)/b],
        [(1/pi)*(0.5-a)*T[12],      T[11]*T[12]/(2*pi**2),   T[12]/(b*pi)],
        [2*(0.5-a),                 T[11]/pi,                2/b]])

    K_non = kap*(U/b)**2*np.array([
        [0., (T[4] + T[10])/pi,               0.],
        [0., (T[5] - T[4]*T[10])/pi**2,       0.],
        [0., 0.,                              0.]])

    K_cir = C*kap*(U/b)**2*np.array([
        [-2*(a+0.5),   -2*(a+0.5)*T[10]/pi,   0.],
        [T[12]/pi,      T[10]*T[12]/pi**2,    0.],
        [2.,            2*T[10]/pi,           0.]])

    return scale @ A_non, scale @ (B_non + B_cir), scale @ (K_non + K_cir)


def pk_flutter(sec, U_vals, structural_damping=True, tol=1e-10, itmax=300):
    """p-k with exact C(k). Returns (U_flutter, f_flutter) or (inf, nan)."""
    b = sec.b
    Cs = sec.Cs_damping() if structural_damping else np.zeros((3, 3))
    w0, _ = sec.modal()
    seeds = [wi*b/U_vals[0] for wi in w0]
    rec = []
    for U in U_vals:
        roots_here = []
        for m in range(3):
            k = seeds[m]
            for _ in range(itmax):
                C = theodorsen(k)
                Ma, Ba, Ka = aero_matrices(sec, U, C)
                M, D, K = sec.Ms + Ma, Cs + Ba, sec.Ks + Ka
                n = 3
                AA = np.block([[np.zeros((n, n)), np.eye(n)],
                               [-inv(M) @ K, -inv(M) @ D]])
                r = eig(AA)[0]
                pos = r[r.imag > 1e-9]
                if len(pos) < 3:
                    k = np.nan; break
                pos = pos[np.argsort(pos.imag)]
                s = pos[m]
                kn = s.imag*b/U
                if abs(kn - k) < tol*max(abs(k), 1.0):
                    k = kn; break
                k = 0.5*k + 0.5*kn
            if np.isnan(k):
                roots_here.append(np.nan + 0j); continue
            C = theodorsen(k)
            Ma, Ba, Ka = aero_matrices(sec, U, C)
            M, D, K = sec.Ms + Ma, Cs + Ba, sec.Ks + Ka
            AA = np.block([[np.zeros((3, 3)), np.eye(3)],
                           [-inv(M) @ K, -inv(M) @ D]])
            r = eig(AA)[0]
            pos = r[r.imag > 1e-9]
            pos = pos[np.argsort(pos.imag)]
            roots_here.append(pos[m])
            seeds[m] = k
        rec.append(roots_here)
    rec = np.array(rec)

    best = (np.inf, np.nan)
    for m in range(3):
        re = rec[:, m].real
        for i in range(len(U_vals)-1):
            if np.isnan(re[i]) or np.isnan(re[i+1]):
                continue
            if re[i] < 0 <= re[i+1]:
                f = -re[i]/(re[i+1]-re[i])
                Uc = U_vals[i] + f*(U_vals[i+1]-U_vals[i])
                wc = abs(rec[i, m]) + f*(abs(rec[i+1, m]) - abs(rec[i, m]))
                if Uc < best[0]:
                    best = (Uc, wc/2/np.pi)
                break
    return best


# ---------------------------------------------------------------------------
# Jones two-lag state space
# ---------------------------------------------------------------------------
# All circulatory aero enters through the single 3/4-chord downwash scalar
#     Q = U*alpha + hdot + b(1/2-a)*alphadot + (T10/pi)*U*beta + (T11/2pi)*b*betadot
# Verified numerically: the circulatory damping and stiffness blocks assembled
# from the TM matrices are both rank 1 and factor exactly as outer(g, Bq) and
# outer(g, Kq) to machine precision, with load direction
#     L = [-(a+1/2), T12/(2pi), 1]
# so the whole thing is g * C(k) * Q. Replacing C(k)*Q by the Wagner/Jones
# indicial form gives the two augmented states Conner refers to.
#
#   Q_eff = 0.5*Q + psi1*eps1*(U/b)*x1 + psi2*eps2*(U/b)*x2
#   xdot_i = -eps_i*(U/b)*x_i + Q
#
# phi(0) = 0.5 is why the instantaneous term carries the 1/2.
#
# STATE ORDERING (8 states, matching Conner's 6 structural + 2 augmented):
#   z = [alpha, beta, h, alphadot, betadot, hdot, x1, x2]

PSI1, EPS1 = 0.165, 0.0455
PSI2, EPS2 = 0.335, 0.3000


def state_space(sec, U, structural_damping=True, freeplay_region="stiff"):
    """A(U) for the 8-state plant. freeplay_region 'stiff' uses the full
    control-surface stiffness; 'free' zeroes K_beta, which is the only thing
    that changes between Conner's piecewise-linear regions."""
    a, b, T = sec.a, sec.b, sec.T
    # noncirculatory pieces: evaluate the assembled matrices at C = 0
    A_non, B_non, K_non = aero_matrices(sec, U, 0.0)
    # circulatory load direction, scaled to dimensional rows
    L = np.array([-(a + 0.5), T[12]/(2*np.pi), 1.0])
    scale = np.diag([sec.M*b**2, sec.M*b**2, sec.M*b])
    g = sec.kappa*(U/b)*(2/b) * (scale @ L)
    Bq = np.array([b*(0.5 - a), T[11]*b/(2*np.pi), 1.0])   # qdot part of Q
    Kq = np.array([U,           T[10]*U/np.pi,     0.0])   # q   part of Q

    Ks = sec.Ks.copy()
    if freeplay_region == "free":
        Ks[1, 1] = 0.0                       # aileron contributes no stiffness
    Cs = sec.Cs_damping() if structural_damping else np.zeros((3, 3))

    Mt = sec.Ms + A_non
    Dt = Cs + B_non + 0.5*np.outer(g, Bq)
    Kt = Ks + K_non + 0.5*np.outer(g, Kq)
    Gl = np.outer(g, [PSI1*EPS1*U/b, PSI2*EPS2*U/b])       # lag feed-through

    Mi = inv(Mt)
    A = np.zeros((8, 8))
    A[0:3, 3:6] = np.eye(3)
    A[3:6, 0:3] = -Mi @ Kt
    A[3:6, 3:6] = -Mi @ Dt
    A[3:6, 6:8] = -Mi @ Gl
    for j, eps in enumerate((EPS1, EPS2)):
        A[6 + j, 0:3] = Kq
        A[6 + j, 3:6] = Bq
        A[6 + j, 6 + j] = -eps*U/b
    return A


def ss_flutter(sec, U_vals, structural_damping=True):
    mx = np.array([np.linalg.eigvals(
        state_space(sec, U, structural_damping)).real.max() for U in U_vals])
    for i in range(len(U_vals) - 1):
        if mx[i] < 0 <= mx[i+1]:
            f = -mx[i]/(mx[i+1] - mx[i])
            Uc = U_vals[i] + f*(U_vals[i+1] - U_vals[i])
            ev = np.linalg.eigvals(state_space(sec, Uc, structural_damping))
            return Uc, abs(ev[np.argmax(ev.real)])/2/np.pi
    return np.inf, np.nan
