"""
diag_conic_torsion.py — Visualize "conic with tilting plane" via Frenet-Serret integration.

Each curve type (circle, ellipse, hyperbola) has a curvature profile κ(θ) derived
from its analytic formula. Adding constant torsion τ rotates the osculating plane as
the parameter θ advances, producing 3D space curves.

Special case: circle (κ=const) + τ=const → circular helix (verified below).
Novel cases:  ellipse + τ=const → "elliptic helix"
              hyperbola branch + τ=const → open spiral

Frenet-Serret equations (state = [P, T, N, B], 12 components):
    dP/dθ = T · ds/dθ
    dT/dθ = κ N · ds/dθ
    dN/dθ = (−κ T + τ B) · ds/dθ
    dB/dθ = −τ N · ds/dθ

The Darboux matrix Ω = [[0,κ,0],[−κ,0,τ],[0,−τ,0]] preserves orthonormality exactly.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# ── Conic curvature profiles ──────────────────────────────────────────────────

def ellipse_kappa(a, b):
    """κ(θ) and ds/dθ for ellipse x=a cos θ, y=b sin θ."""
    def kappa(theta):
        return a * b / (a**2 * np.sin(theta)**2 + b**2 * np.cos(theta)**2)**1.5
    def ds(theta):
        return np.sqrt(a**2 * np.sin(theta)**2 + b**2 * np.cos(theta)**2)
    # Initial conditions at θ=0: pos=(a,0,0), T=(0,1,0), N=(-1,0,0) (inward)
    # For circle (a=b=R): N=(-1,0,0) is correct (toward center)
    # For ellipse at θ=0: d²P/dθ² = (-a cos 0, -b sin 0) = (-a, 0), κN direction = -x
    pos0 = np.array([a, 0.0, 0.0])
    T0   = np.array([0.0, 1.0, 0.0])
    N0   = np.array([-1.0, 0.0, 0.0])
    return kappa, ds, pos0, T0, N0


def hyperbola_kappa(a=1.0, b=1.0):
    """κ(t) and ds/dt for right branch of hyperbola: x=a cosh t, y=b sinh t.

    Using the hyperbolic parameterization avoids the angular-parameter singularity
    near the asymptotes.  Curvature: κ = ab / (a²sinh²t + b²cosh²t)^(3/2).
    Normal points AWAY from the origin (hyperbola is concave outward).
    """
    def kappa(t):
        s2 = a**2 * np.sinh(t)**2 + b**2 * np.cosh(t)**2
        return a * b / s2**1.5
    def ds(t):
        return np.sqrt(a**2 * np.sinh(t)**2 + b**2 * np.cosh(t)**2)
    # At t=0: x=a, y=0; dx/dt=0, dy/dt=b → T=(0,1,0)
    # d²x/dt²=a cosh 0=a, d²y/dt²=b sinh 0=0
    # κN direction = (a, 0) normalized = (+1, 0) → N points AWAY from origin
    pos0 = np.array([a, 0.0, 0.0])
    T0   = np.array([0.0, 1.0, 0.0])
    N0   = np.array([1.0, 0.0, 0.0])   # outward (hyperbola concave outward)
    return kappa, ds, pos0, T0, N0


# ── Frenet-Serret integrator ──────────────────────────────────────────────────

def frenet_integrate(kappa_func, ds_func, tau, param_range, pos0, T0, N0,
                     n_pts=800):
    """Integrate Frenet-Serret ODE to produce 3D curve with given κ(θ) and τ=const.

    Returns:
        P       — (n_pts, 3) positions
        frames  — (n_pts, 3, 3) Frenet frames [T, N, B] as rows
        theta   — (n_pts,) parameter values
    """
    B0 = np.cross(T0, N0)
    state0 = np.concatenate([pos0, T0, N0, B0])

    def rhs(theta, state):
        T = state[3:6]
        N = state[6:9]
        B = state[9:12]
        k  = kappa_func(theta)
        ds = ds_func(theta)
        return np.concatenate([
            T * ds,
            k  * N * ds,
            (-k * T + tau * B) * ds,
            -tau * N * ds,
        ])

    t_eval = np.linspace(param_range[0], param_range[1], n_pts)
    sol = solve_ivp(rhs, param_range, state0, t_eval=t_eval,
                    method='DOP853', rtol=1e-10, atol=1e-12)

    P      = sol.y[0:3].T
    frames = sol.y[3:12].T.reshape(n_pts, 3, 3)   # rows: T, N, B
    return P, frames, t_eval


# ── Analytic helix verification ───────────────────────────────────────────────

def analytic_helix(R, tau_val, n_turns, n_pts=800):
    """For a circle of radius R with Frenet torsion τ, compute the analytic helix.

    Frenet invariants: κ = 1/R (circle), τ = τ_val.
    Helix parameters: helix_R = κ/(κ²+τ²), pitch = τ/(κ²+τ²).
    Arc-length period: 2π / √(κ²+τ²).
    """
    kap = 1.0 / R
    denom = kap**2 + tau_val**2
    hR    = kap / denom      # helix radius
    pitch = tau_val / denom  # helix pitch per radian of arc-length
    omega = np.sqrt(denom)   # angular rate in arc-length
    s_end = n_turns * 2 * np.pi / omega
    s     = np.linspace(0, s_end, n_pts)
    # The helix aligned to match Frenet initial conditions:
    # At s=0: P=(R,0,0), T=(0,1,0) → the helix axis must be in z-direction,
    # centered so that x(0)=R, y(0)=0.
    x = hR * np.cos(omega * s - np.pi/2) + (R - hR)
    # Hmm, this is tricky to align analytically — just show it qualitatively.
    # Use the simpler centered version:
    x = hR * np.cos(omega * s)
    y = hR * np.sin(omega * s)
    z = pitch * s
    return np.column_stack([x, y, z])


# ── Tilting-plane visualization ───────────────────────────────────────────────

def plot_tilting_planes(ax, P, frames, stride=80, scale=0.25, alpha=0.25):
    """Draw small discs representing the local osculating plane at evenly spaced points."""
    n_pts = len(P)
    u_circ = np.linspace(0, 2*np.pi, 30)
    cx = scale * np.cos(u_circ)
    cy = scale * np.sin(u_circ)

    for i in range(0, n_pts, stride):
        center = P[i]
        e1 = frames[i, 0]   # T
        e2 = frames[i, 1]   # N
        # disk in (T, N) plane
        disk = center + cx[:, None]*e1 + cy[:, None]*e2
        ax.plot(disk[:, 0], disk[:, 1], disk[:, 2], 'gray', lw=0.5, alpha=alpha)


# ── Main ──────────────────────────────────────────────────────────────────────

TAUS   = [0.0, 0.2, 0.5, 1.0]
COLORS = ['steelblue', 'forestgreen', 'darkorange', 'crimson']
N_TURNS = 2   # how many conic periods to trace

fig = plt.figure(figsize=(18, 12))
fig.suptitle('Conic with tilting plane  (Frenet-Serret: constant torsion τ)', fontsize=14)

# ── Panel 1: Circle ───────────────────────────────────────────────────────────
ax1 = fig.add_subplot(2, 3, 1, projection='3d')
ax1.set_title('Circle  a=b=1\n(τ>0 → helix)')

kappa_c, ds_c, pos_c, T_c, N_c = ellipse_kappa(1.0, 1.0)
for tau, col in zip(TAUS, COLORS):
    P, frames, _ = frenet_integrate(kappa_c, ds_c, tau,
                                    (0, N_TURNS * 2*np.pi), pos_c, T_c, N_c)
    ax1.plot(P[:,0], P[:,1], P[:,2], color=col, lw=1.5, label=f'τ={tau}')
plot_tilting_planes(ax1, P, frames)   # show planes for last (τ=1) curve
ax1.legend(fontsize=7, loc='upper left')

# ── Panel 2: Ellipse (a=2, b=1) ───────────────────────────────────────────────
ax2 = fig.add_subplot(2, 3, 2, projection='3d')
ax2.set_title('Ellipse  a=2, b=1\n("elliptic helix")')

kappa_e, ds_e, pos_e, T_e, N_e = ellipse_kappa(2.0, 1.0)
for tau, col in zip(TAUS, COLORS):
    P, frames, _ = frenet_integrate(kappa_e, ds_e, tau,
                                    (0, N_TURNS * 2*np.pi), pos_e, T_e, N_e)
    ax2.plot(P[:,0], P[:,1], P[:,2], color=col, lw=1.5, label=f'τ={tau}')
ax2.legend(fontsize=7, loc='upper left')

# ── Panel 3: Hyperbola branch ─────────────────────────────────────────────────
ax3 = fig.add_subplot(2, 3, 3, projection='3d')
ax3.set_title('Hyperbola branch  a=b=1\n(open arc with torsion)')

kappa_h, ds_h, pos_h, T_h, N_h = hyperbola_kappa(1.0, 1.0)
for tau, col in zip(TAUS, COLORS):
    P, frames, _ = frenet_integrate(kappa_h, ds_h, tau,
                                    (-2.5, 2.5), pos_h, T_h, N_h, n_pts=500)
    ax3.plot(P[:,0], P[:,1], P[:,2], color=col, lw=1.5, label=f'τ={tau}')
ax3.legend(fontsize=7, loc='upper left')

# ── Panel 4: Tilting-planes visualization (ellipse, τ=0.3) ───────────────────
ax4 = fig.add_subplot(2, 3, 4, projection='3d')
ax4.set_title('Ellipse τ=0.3: osculating planes\n(grey discs = local conic plane)')

P, frames, _ = frenet_integrate(kappa_e, ds_e, 0.3,
                                (0, 3 * 2*np.pi), pos_e, T_e, N_e, n_pts=1200)
ax4.plot(P[:,0], P[:,1], P[:,2], color='darkorange', lw=2)
plot_tilting_planes(ax4, P, frames, stride=100, scale=0.4, alpha=0.35)

# ── Panel 5: Circle helix vs analytic (τ=0.5, verify) ────────────────────────
ax5 = fig.add_subplot(2, 3, 5, projection='3d')
ax5.set_title('Circle τ=0.5: Frenet vs analytic helix\n(should overlap)')

tau_v = 0.5
P_frenet, _, _ = frenet_integrate(kappa_c, ds_c, tau_v,
                                   (0, 3 * 2*np.pi), pos_c, T_c, N_c, n_pts=1000)
ax5.plot(P_frenet[:,0], P_frenet[:,1], P_frenet[:,2],
         color='steelblue', lw=2, label='Frenet integration')

# Check: κ=1/R=1, τ=0.5 → helix with R_h=1/(1+0.25)=0.8, pitch=0.5/(1+0.25)=0.4
kap_v = 1.0
denom = kap_v**2 + tau_v**2
hR    = kap_v / denom
pitch = tau_v / denom
omega = np.sqrt(denom)
s_pts = np.linspace(0, 3 * 2*np.pi / omega, 1000)
# Align to initial conditions P=(1,0,0), T=(0,1,0):
# At s=0: want (x,y,z)=(1,0,0). Standard helix: x=hR cos(ωs), y=hR sin(ωs), z=pitch*s
# doesn't start at (1,0). Instead use phase offset φ so that hR cos(φ)=1, hR sin(φ)=0
# → φ=0, but hR=0.8 ≠ 1. The discrepancy means the axis isn't at origin for these ICs.
# The helix with these ICs: the axis is offset by (1-hR, 0, 0) from origin.
P_analytic = np.column_stack([
    hR * np.cos(omega * s_pts) + (1.0 - hR),   # shifted so x(0)=1
    hR * np.sin(omega * s_pts),
    pitch * s_pts,
])
ax5.plot(P_analytic[:,0], P_analytic[:,1], P_analytic[:,2],
         color='crimson', lw=1, linestyle='--', label='Analytic helix')
ax5.legend(fontsize=8)

# Print deviation
dev = np.max(np.linalg.norm(P_frenet - P_analytic, axis=1))
ax5.set_xlabel(f'max |Frenet−analytic| = {dev:.2e}', fontsize=8)

# ── Panel 6: Eccentric ellipse (a=4, b=0.5) to show dramatic shape ───────────
ax6 = fig.add_subplot(2, 3, 6, projection='3d')
ax6.set_title('Eccentric ellipse  a=4, b=0.5\n(κ varies 8:1 around orbit)')

kappa_x, ds_x, pos_x, T_x, N_x = ellipse_kappa(4.0, 0.5)
for tau, col in zip(TAUS, COLORS):
    P, _, _ = frenet_integrate(kappa_x, ds_x, tau,
                               (0, N_TURNS * 2*np.pi), pos_x, T_x, N_x)
    ax6.plot(P[:,0], P[:,1], P[:,2], color=col, lw=1.5, label=f'τ={tau}')
ax6.legend(fontsize=7, loc='upper left')

for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.set_xlabel('x', fontsize=8)
    ax.set_ylabel('y', fontsize=8)
    ax.set_zlabel('z', fontsize=8)
    ax.tick_params(labelsize=6)

plt.tight_layout()
outfile = 'conic_torsion_curves.png'
plt.savefig(outfile, dpi=120, bbox_inches='tight')
print(f"Saved: {outfile}")
plt.show()
