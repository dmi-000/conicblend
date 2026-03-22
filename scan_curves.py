#!/usr/bin/env python3
"""scan_curves.py — visualise demo_scan output.

Compiles and runs demo_scan if CSV files are missing, then produces:
  demo_blend.png   — 2×4 grid: curve view (top) + convergence (bottom)
"""

import os
import subprocess
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.collections as mc

# ── compile / run demo_scan ────────────────────────────────────────────────

BINARY  = os.path.join(os.path.dirname(__file__), "demo_scan")
SRC     = BINARY + ".cpp"
SDK_INC = "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1"
INC_DIR = os.path.dirname(__file__)

def _need_build():
    if not os.path.exists(BINARY):
        return True
    return os.path.getmtime(SRC) > os.path.getmtime(BINARY)

def _build():
    flags = ["-std=c++17", "-O2", f"-I{INC_DIR}"]
    if os.path.isdir(SDK_INC):
        flags.append(f"-I{SDK_INC}")
    cmd = ["g++", *flags, "-o", BINARY, SRC]
    print("Compiling demo_scan …", " ".join(cmd))
    subprocess.check_call(cmd)

def ensure_csvs():
    if _need_build():
        _build()
    # Run in frenetconic/ so CSVs land there
    cwd = os.path.dirname(os.path.abspath(__file__))
    missing = not os.path.exists(os.path.join(cwd, "ellipse_conv.csv"))
    if missing:
        print("Running demo_scan …")
        subprocess.check_call([BINARY], cwd=cwd)

# ── CSV loaders ────────────────────────────────────────────────────────────

DIR = os.path.dirname(os.path.abspath(__file__))

def load(name):
    return np.genfromtxt(os.path.join(DIR, name), delimiter=",", names=True)

# ── plot helpers ───────────────────────────────────────────────────────────

C_TRUTH  = "#aaaaaa"
C_CONIC  = "#1f77b4"   # blue
C_CIRCLE = "#ff7f0e"   # orange
C_CTRL   = "#333333"
C_VALID  = "#2ca02c"   # green  (valid window centre)
C_INVAL  = "#d62728"   # red    (invalid window centre)

def plot_curve(ax, name, title):
    truth  = load(f"{name}_truth.csv")
    blend  = load(f"{name}_blend.csv")
    circle = load(f"{name}_circle.csv")
    ctrl   = load(f"{name}_ctrl.csv")
    wins   = load(f"{name}_wins.csv")

    ax.plot(truth["x"],  truth["y"],  color=C_TRUTH,  lw=2.5, zorder=1, label="truth")
    ax.plot(circle["x"], circle["y"], color=C_CIRCLE, lw=1.5, ls="--",
            zorder=2, label="circle")
    ax.plot(blend["x"],  blend["y"],  color=C_CONIC,  lw=1.5, zorder=3, label="conic")
    ax.scatter(ctrl["x"], ctrl["y"],  s=18, color=C_CTRL, zorder=5, label="ctrl")

    # Window centre marks: green = valid, red = invalid
    for row in wins:
        if np.isscalar(row):  # single-row edge case
            break
        # wins array has fields i, t_mid, valid
    # Use structured array directly
    valid_mask = wins["valid"].astype(bool)
    # Need x,y positions of t_mid from ctrl (t_mid ≈ ctrl[i+2].t)
    ctrl_t = ctrl["t"]
    ctrl_x = ctrl["x"]
    ctrl_y = ctrl["y"]
    for row in wins:
        i_win = int(row["i"])
        # Window i covers ctrl[i..i+4]; centre is ctrl[i+2]
        idx = i_win + 2
        if idx < len(ctrl_t):
            col = C_VALID if row["valid"] else C_INVAL
            ax.scatter(ctrl_x[idx], ctrl_y[idx], s=50, color=col,
                       zorder=6, marker="^", linewidths=0)

    ax.set_title(title, fontsize=10)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=7)

def plot_convergence(ax, name, title):
    data = load(f"{name}_conv.csv")
    ns     = data["n"]
    ec     = data["err_conic"]
    ek     = data["err_circle"]
    pct    = data["pct_valid_wins"]

    ax.loglog(ns, ek, color=C_CIRCLE, lw=1.4, ls="--", marker="o",
              ms=4, label="circle")
    ax.loglog(ns, ec, color=C_CONIC,  lw=1.4, ls="-",  marker="s",
              ms=4, label="conic")

    # h^2 reference line through circle's last point
    n_ref = np.array([ns[0], ns[-1]], dtype=float)
    scale = ek[-1] * (ns[-1] ** 2)
    ax.loglog(n_ref, scale / n_ref**2, color="gray", lw=0.8, ls=":", label="$n^{-2}$")

    ax.set_xlabel("n", fontsize=8)
    ax.set_ylabel("max error", fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, loc="upper right")
    ax.tick_params(labelsize=7)

    # annotate max % valid windows
    ax.text(0.03, 0.06,
            f"valid wins: {pct.min():.0f}–{pct.max():.0f}%",
            transform=ax.transAxes, fontsize=7, color="#555555")

# ── main ───────────────────────────────────────────────────────────────────

CURVES = [
    ("ellipse",   "Tilted ellipse"),
    ("parabola",  "Parabola  y = x²"),
    ("lissajous", "Lissajous 3:2"),
    ("rose",      "Rose  k = 3"),
]

def main():
    ensure_csvs()

    fig, axes = plt.subplots(2, 4, figsize=(14, 7))
    fig.suptitle(
        "conicblend demo  —  conic (blue) vs circle (orange dashed) blend\n"
        "triangles: ▲ valid window  ▲ invalid window  (green / red)",
        fontsize=10)

    for col, (name, title) in enumerate(CURVES):
        plot_curve      (axes[0][col], name, title)
        plot_convergence(axes[1][col], name, title + " — convergence")

    # Shared legend for row 0
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles[:4], labels[:4],
               loc="upper right", fontsize=8, ncol=2,
               bbox_to_anchor=(0.99, 0.97))

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(DIR, "demo_blend.png")
    plt.savefig(out, dpi=130, bbox_inches="tight")
    print(f"Saved: {out}")

if __name__ == "__main__":
    main()
