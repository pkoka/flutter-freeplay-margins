"""
Direct port of Hodges' Maple worksheet CA2.mws / H5_dec3.mws.

Independent of typical_section.py: his aero enters as nondimensional
coefficients lh, ltheta, mh, mtheta from the Bisplinghoff-style tables,
not as the dimensional Theodorsen L and M I wrote. Different algebra,
same physics. If the flutter points agree, my derivation is clean.

His eigenvalue is  x = (w_theta/w)^2 (1 + i g),  so
  w/w_theta = 1/sqrt(Re x),  g = Im x / Re x,  U/(b w_theta) = (w/w_theta)/k
"""
import numpy as np
from scipy.special import hankel2

def C_exact(k):
    H1, H0 = hankel2(1, k), hankel2(0, k)
    return H1/(H1 + 1j*H0)

def C_rational(k):          # the approximation quoted in H&P Example 2
    return (0.01365 + 0.2808j*k - k**2/2)/(0.01365 + 0.3455j*k - k**2)

def sweep(a, e, mu, r2, sigma, k_vals, Cf=C_exact):
    Xth = e - a
    out = []
    for k in k_vals:
        C = Cf(k)
        lh     = 1 - 2j*C/k
        lth    = -a - 1j/k - 2*C/k**2 - 2j*(0.5 - a)*C/k
        mh     = -a + 2j*(0.5 + a)*C/k
        mth    = 1/8 + a**2 - 1j*(0.5 - a)/k + 2*(0.5 + a)*C/k**2 \
                 + 2j*(0.25 - a**2)*C/k
        # det A = 0 is quadratic in x
        A11_0, A11_1 = mu + lh,        -mu*sigma**2
        A22_0, A22_1 = mu*r2 + mth,    -mu*r2
        B12, B21 = mu*Xth + lth,       mu*Xth + mh
        qa = A11_1*A22_1
        qb = A11_0*A22_1 + A11_1*A22_0
        qc = A11_0*A22_0 - B12*B21
        for x in np.roots([qa, qb, qc]):
            if x.real > 1e-12:
                w = 1/np.sqrt(x.real)
                out.append((k, w/k, w, x.imag/x.real))   # k, V*, w/w_th, g
    return np.array(out)

def flutter(rows):
    """lowest V* where g crosses zero from below, per branch"""
    best = (np.inf, np.nan)
    for w0, w1 in ((0.0, 0.75), (0.75, 3.0)):            # crude branch split
        s = rows[(rows[:, 2] > w0) & (rows[:, 2] <= w1)]
        s = s[np.argsort(s[:, 1])]
        for i in range(len(s)-1):
            if s[i, 3] < 0 <= s[i+1, 3]:
                f = -s[i, 3]/(s[i+1, 3] - s[i, 3])
                V = s[i, 1] + f*(s[i+1, 1] - s[i, 1])
                W = s[i, 2] + f*(s[i+1, 2] - s[i, 2])
                if V < best[0]:
                    best = (V, W)
                break
    return best

k = np.geomspace(2.0, 0.02, 20000)

print("CA2.mws  Example 1:  a=-1/5  e=-1/10  mu=20  r2=6/25  sigma=2/5")
V, W = flutter(sweep(-1/5, -1/10, 20, 6/25, 2/5, k))
print(f"   Hodges k-method, exact C(k)   V* = {V:.4f}   w_f/w_th = {W:.4f}")

print("\nH5_dec3.mws  Example 2:  a=-1/10  e=-1/15  mu=25  r2=7/25  sigma=2/5")
V, W = flutter(sweep(-1/10, -1/15, 25, 7/25, 2/5, k))
print(f"   exact C(k)                    V* = {V:.4f}   w_f/w_th = {W:.4f}")
V, W = flutter(sweep(-1/10, -1/15, 25, 7/25, 2/5, k, C_rational))
print(f"   rational C(k) approximation   V* = {V:.4f}   w_f/w_th = {W:.4f}")
