"""
fit_frenet_5pt.py — Fit κ and τ from 5 3D sample points.

Algorithm (one 5-point window):
  1. SVD best-fit plane → local frame (e1≈T, e2≈N, e3=e1×e2≈B)
  2. Menger curvature κ from the 3 inner projected points
  3. Cumulative chord lengths → arc-length positions s₀…s₄ (s₂=0)
  4. Least-squares  w(s) ≈ A·s³  using the 5 out-of-plane residuals
  5. τ = 6A / κ

Why s³? Frenet expansion:
    P(s) = P₀ + s·T + (κs²/2)·N + (κτs³/6)·B + O(s⁴)
    ⟹ w = (P-P₀)·B = (κτ/6)s³   (leading B-component is cubic in arc-length)
"""

import numpy as np
import matplotlib.pyplot as plt


# ── Core estimator ────────────────────────────────────────────────────────────

def fit_kappa_tau(pts5):
    """Estimate (κ, τ) from 5 3D points.

    Uses the Frenet frame at pts5[2] computed from discrete differences —
    NOT the SVD plane, which fails when torsion is large (helix z-drift
    becomes a principal SVD direction, corrupting the w-residuals).

    Frenet frame at center (pts5[2]):
        T = central-difference tangent  (P3 - P1) / |P3 - P1|
        N = normal component of second difference  (P1 - 2P2 + P3)
        B = T × N

    Returns
    -------
    kappa : float  — Menger curvature at center point pts5[2]
    tau   : float  — estimated torsion
    w     : (5,)   — B-component of (pts5 - pts5[2])
    s     : (5,)   — arc-length positions, s[2] = 0
    A     : float  — fitted cubic coefficient (= κτ/6)
    """
    p = pts5

    # 1. Frenet frame at p[2] via central differences
    T_raw = p[3] - p[1]
    T_norm = np.linalg.norm(T_raw)
    if T_norm < 1e-15:
        return 0.0, 0.0, np.zeros(5), np.zeros(5), 0.0
    T = T_raw / T_norm

    acc = p[1] - 2*p[2] + p[3]        # ∝ d²P/dθ² · Δθ²  (acceleration)
    acc_N = acc - np.dot(acc, T) * T   # remove tangential part → normal direction
    N_norm = np.linalg.norm(acc_N)
    if N_norm < 1e-15:
        return 0.0, 0.0, np.zeros(5), np.zeros(5), 0.0
    N = acc_N / N_norm
    B = np.cross(T, N)                 # binormal: (T, N, B) right-handed

    # 2. Project all 5 points into Frenet frame at p[2]
    ctr = p - p[2]
    w   = ctr @ B                      # B-component: should follow (κτ/6)s³

    # 3. Arc-length positions (cumulative chord lengths, centered at index 2)
    chords = np.linalg.norm(np.diff(p, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(chords)])
    s -= s[2]

    # 4. Menger curvature at p[2] from inner 3 points
    a = np.linalg.norm(p[2] - p[1])
    b = np.linalg.norm(p[3] - p[2])
    c = np.linalg.norm(p[3] - p[1])
    area  = 0.5 * np.linalg.norm(np.cross(p[2] - p[1], p[3] - p[1]))
    kappa = 2 * area / (a * b * c) if a * b * c > 1e-15 else 0.0

    # 5. Fit  w ≈ A·s³  (1-parameter least squares: A = Σwᵢsᵢ³ / Σsᵢ⁶)
    s3 = s**3
    A  = np.dot(w, s3) / (np.dot(s3, s3) + 1e-30)

    tau = 6.0 * A / kappa if abs(kappa) > 1e-10 else 0.0
    return kappa, tau, w, s, A


# ── Test curves ───────────────────────────────────────────────────────────────

def helix_sample(R, pitch, theta_center, delta_theta):
    """5 points on helix x=R cosθ, y=R sinθ, z=pitch·θ.

    True Frenet invariants:
        L = sqrt(R²+p²),  κ = R/L²,  τ = p/L²
    """
    L           = np.sqrt(R**2 + pitch**2)
    kappa_true  = R / L**2
    tau_true    = pitch / L**2
    thetas      = theta_center + delta_theta * np.array([-2, -1, 0, 1, 2])
    pts5        = np.column_stack([R*np.cos(thetas), R*np.sin(thetas), pitch*thetas])
    return pts5, kappa_true, tau_true


def torus_knot_sample(t_center, delta_t, p=2, q=3, R=2.0, r=1.0):
    """5 points on a (p,q) torus knot."""
    ts    = t_center + delta_t * np.array([-2, -1, 0, 1, 2])
    phi   = p * ts;  theta = q * ts
    pts5  = np.column_stack([
        (R + r*np.cos(theta)) * np.cos(phi),
        (R + r*np.cos(theta)) * np.sin(phi),
        r * np.sin(theta),
    ])
    return pts5


# ── Main ──────────────────────────────────────────────────────────────────────

R, pitch = 1.0, 0.5
L_h = np.sqrt(R**2 + pitch**2)
kappa_true = R / L_h**2
tau_true   = pitch / L_h**2

print(f"Helix R={R}, pitch={pitch}")
print(f"  True: κ = {kappa_true:.6f},  τ = {tau_true:.6f}")

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle(
    r'Fitting $\kappa$ and $\tau$ from 5 3D points'
    '\n'
    r'Frenet expansion: $w(s) \approx \frac{\kappa\tau}{6}\,s^3$'
    r'  $\Rightarrow$  $\tau = 6A/\kappa$',
    fontsize=13)


# ── Panel 1: w vs s for single helix window ───────────────────────────────────
ax = axes[0, 0]
pts5, _, _ = helix_sample(R, pitch, np.pi/4, 0.35)
k_rec, t_rec, w_coords, s_coords, A_fit = fit_kappa_tau(pts5)
print(f"  Δθ=0.35:  κ_rec={k_rec:.6f}, τ_rec={t_rec:.6f}")

s_dense = np.linspace(s_coords[0]*1.1, s_coords[-1]*1.1, 100)
ax.plot(s_dense, A_fit * s_dense**3, 'r--', lw=1.5,
        label=rf'fit $A\cdot s^3$, $\tau_{{rec}}$={t_rec:.4f}')
ax.scatter(s_coords, w_coords, color='steelblue', s=70, zorder=5,
           label='w data')
ax.axhline(0, color='gray', lw=0.5)
ax.axvline(0, color='gray', lw=0.5)
ax.set_xlabel('arc-length s')
ax.set_ylabel('w  (out-of-plane)')
ax.set_title(f'Helix R={R}, p={pitch}\n'
             rf'True $\tau$={tau_true:.4f},  Rec $\tau$={t_rec:.4f}')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)


# ── Panel 2: Convergence of κ, τ vs window size ───────────────────────────────
ax = axes[0, 1]
deltas = np.geomspace(0.01, 1.8, 40)
kappas_rec, taus_rec = [], []
for d in deltas:
    pts5d, _, _ = helix_sample(R, pitch, np.pi/4, d)
    k_d, t_d, *_ = fit_kappa_tau(pts5d)
    kappas_rec.append(k_d)
    taus_rec.append(t_d)
kappas_rec = np.array(kappas_rec)
taus_rec   = np.array(taus_rec)

ax.loglog(deltas, np.abs(taus_rec   - tau_true)   / tau_true   + 1e-16,
          'steelblue', lw=2, label=r'$|\tau_{rec}-\tau_{true}|/\tau_{true}$')
ax.loglog(deltas, np.abs(kappas_rec - kappa_true) / kappa_true + 1e-16,
          'darkorange', lw=2, label=r'$|\kappa_{rec}-\kappa_{true}|/\kappa_{true}$')
# Reference slopes
ax.loglog(deltas, 0.3 * deltas**2, 'gray', lw=1, ls='--', label=r'$O(\Delta\theta^2)$')
ax.loglog(deltas, 0.05 * deltas**4, 'gray', lw=1, ls=':',  label=r'$O(\Delta\theta^4)$')
ax.set_xlabel(r'$\Delta\theta$ (window half-spacing, rad)')
ax.set_ylabel('Relative error')
ax.set_title(r'Convergence of $\kappa$, $\tau$ estimates')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)


# ── Panel 3: τ recovery across a range of pitches ────────────────────────────
ax = axes[0, 2]
pitches_v = np.linspace(0.0, 2.5, 50)
taus_true_v, taus_rec_v = [], []
for pv in pitches_v:
    Lv = np.sqrt(1.0 + pv**2)
    taus_true_v.append(pv / Lv**2)
    pts5v, _, _ = helix_sample(1.0, pv, np.pi/4, 0.3)
    _, t_v, *_ = fit_kappa_tau(pts5v)
    taus_rec_v.append(t_v)

ax.plot(pitches_v, taus_true_v, 'k-',  lw=2,  label=r'True $\tau$')
ax.scatter(pitches_v, taus_rec_v, color='steelblue', s=18, zorder=5,
           label=r'Recovered $\tau$  ($\Delta\theta$=0.3)')
ax.set_xlabel('Helix pitch p')
ax.set_ylabel(r'$\tau$')
ax.set_title(r'$\tau$ recovery across pitch range''\n(κ varies too: max τ at p≈R)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)


# ── Panel 4: Torus knot — κ(t) and τ(t) around the orbit ────────────────────
ax = axes[1, 0]
t_centers = np.linspace(0.15, 2*np.pi - 0.15, 60)
kaps_k, taus_k = [], []
for tc in t_centers:
    pts5k = torus_knot_sample(tc, 0.08)
    k_k, t_k, *_ = fit_kappa_tau(pts5k)
    kaps_k.append(k_k)
    taus_k.append(t_k)

ax2 = ax.twinx()
ln1, = ax.plot(t_centers, kaps_k, 'darkorange', lw=1.8, label=r'$\kappa$ (left)')
ln2, = ax2.plot(t_centers, taus_k, 'steelblue',  lw=1.8, label=r'$\tau$ (right)')
ax.set_xlabel('t  (parameter around torus knot)')
ax.set_ylabel(r'Menger $\kappa$', color='darkorange')
ax2.set_ylabel(r'Estimated $\tau$', color='steelblue')
ax.set_title('(2,3) Torus knot: recovered κ(t), τ(t)\n'
             r'($\Delta t$=0.08, R=2, r=1)')
ax.legend(handles=[ln1, ln2], fontsize=8)
ax.grid(True, alpha=0.3)


# ── Panel 5: Planar curve — τ should vanish ───────────────────────────────────
ax = axes[1, 1]
# Eccentric planar ellipse (a=3, b=1, z=0)
thetas_e = np.linspace(0, 2*np.pi, 200, endpoint=False)
taus_planar = []
thetas_mid  = []
for i in range(2, len(thetas_e) - 2):
    idxs  = [i-2, i-1, i, i+1, i+2]
    pts5p = np.column_stack([
        3 * np.cos(thetas_e[idxs]),
        1 * np.sin(thetas_e[idxs]),
        np.zeros(5),
    ])
    _, t_p, *_ = fit_kappa_tau(pts5p)
    taus_planar.append(t_p)
    thetas_mid.append(thetas_e[i])

taus_planar = np.array(taus_planar)
noise = np.std(taus_planar)
ax.plot(thetas_mid, taus_planar, 'steelblue', lw=1.2)
ax.axhline(0, color='k', lw=1, ls='--', label=r'True $\tau$ = 0')
ax.set_xlabel(r'$\theta$ around ellipse')
ax.set_ylabel(r'Recovered $\tau$  (should be $\approx 0$)')
ax.set_title(f'Planar ellipse (a=3, b=1):  τ noise floor\n'
             rf'std($\tau_{{rec}}$) = {noise:.1e}  (machine precision)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)


# ── Panel 6: w fit quality at 3 window sizes ─────────────────────────────────
ax = axes[1, 2]
COLORS3 = ['steelblue', 'darkorange', 'crimson']
for dth, col in zip([0.1, 0.4, 1.0], COLORS3):
    pts5h, _, _ = helix_sample(R, pitch, np.pi/4, dth)
    k_h, t_h, w_h, s_h, A_h = fit_kappa_tau(pts5h)
    err_pct = 100 * abs(t_h - tau_true) / tau_true
    s_d = np.linspace(s_h[0]*1.05, s_h[-1]*1.05, 80)
    ax.plot(s_d, A_h * s_d**3, '--', color=col, lw=1.2)
    ax.scatter(s_h, w_h, color=col, s=55, zorder=5,
               label=rf'$\Delta\theta$={dth}: $\tau$={t_h:.3f} ({err_pct:.1f}% err)')

ax.axhline(0, color='gray', lw=0.5)
ax.axvline(0, color='gray', lw=0.5)
ax.set_xlabel('arc-length s')
ax.set_ylabel('w  (out-of-plane)')
ax.set_title(rf'$w \approx A\cdot s^3$ quality at 3 window sizes'
             f'\nTrue τ={tau_true:.4f}')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)


plt.tight_layout()
outfile = 'fit_frenet_5pt.png'
plt.savefig(outfile, dpi=120, bbox_inches='tight')
print(f"\nSaved: {outfile}")
plt.show()
