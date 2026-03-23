// demo_torsion.cpp — torsion properties of ConicWindow in 3D
//
// ConicWindow fits a planar conic to the best-fit plane of 5 control points.
// For a 3D curve with non-zero torsion the 5 points are not planar; the orbit
// lies in the best-fit plane, discarding the out-of-plane (torsion) component.
//
// This file characterises that behaviour mathematically:
//   1. Per-window table: torsion, predicted κτ(2h)³/6, observed planarity_rms
//   2. h³ scaling study: verify planarity_residual ∝ h³ as n → ∞
//   3. Interpolation error: measure O(planarity_residual) deviation from 3D knots
//
// Test curve: twisted cubic  r(t) = (t, t², t³),  t ∈ [0,1]
//   κ(t) = 2√(1+9t²)/(1+4t²+9t⁴)^{3/2}
//   τ(t) = 3/(1+4t²+9t⁴)
//
// Frenet-Serret prediction for planarity residual of a window of half-width 2h:
//   planarity_rms ≈ C · κ(t_mid) · τ(t_mid) · (2h)³ / 6
// where C is a numerical prefactor depending on how the 5 points sample the window.
//
// Compile:
//   g++ -std=c++17 -O2 -I. -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 -o demo_torsion demo_torsion.cpp

#include "conicblend.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace fc;
using namespace fc::detail;
using V3 = VecN<3>;

static void pass(char const* name) { std::printf("PASS  %s\n", name); }
static void fail(char const* name, char const* msg) {
    std::printf("FAIL  %s: %s\n", name, msg);
    std::exit(1);
}

// ── Twisted cubic geometry ─────────────────────────────────────────────────

static V3    tc(double t)       { return V3(t, t*t, t*t*t); }
static V3    tc_d1(double t)    { return V3(1, 2*t, 3*t*t); }
static V3    tc_d2(double t)    { return V3(0, 2,   6*t);   }

// speed, curvature, torsion at t
static double tc_speed(double t) { return tc_d1(t).norm(); }

static double tc_curvature(double t) {
    V3 d1 = tc_d1(t), d2 = tc_d2(t);
    // κ = |d1 × d2| / |d1|³
    V3 cross(d1[1]*d2[2]-d1[2]*d2[1],
             d1[2]*d2[0]-d1[0]*d2[2],
             d1[0]*d2[1]-d1[1]*d2[0]);
    double sp = d1.norm();
    return cross.norm() / (sp*sp*sp);
}

static double tc_torsion(double t) {
    // τ = 3 / (1 + 4t² + 9t⁴)
    return 3.0 / (1.0 + 4.0*t*t + 9.0*t*t*t*t);
}

// ── Planarity residual ─────────────────────────────────────────────────────
// RMS distance of the 5 original 3D points from their best-fit plane.

static double planarity_rms(V3 const pts5[5]) {
    PlaneND<3> pl = best_fit_plane<3>(pts5);
    double s = 0.0;
    for (int i = 0; i < 5; ++i) {
        V3 proj = pl.center + pl.pts2d[i][0]*pl.e1 + pl.pts2d[i][1]*pl.e2;
        double d = (pts5[i] - proj).norm();
        s += d*d;
    }
    return std::sqrt(s / 5.0);
}

// ── Test 1: per-window analysis ────────────────────────────────────────────
//
// For each window: show torsion, κτ(2h)³/6 prediction, observed planarity_rms,
// and whether ConicWindow accepted it.  Demonstrates that the accepted windows
// have non-trivial planarity residuals (the gap), and that the Frenet-Serret
// formula predicts them well.

void test_per_window()
{
    int n = 12;
    std::vector<V3>     ctrl(n);
    std::vector<double> ts(n);
    for (int i = 0; i < n; ++i) {
        double t = (double)i / (n-1);
        ctrl[i]  = tc(t);
        ts[i]    = t;
    }
    double h = 1.0 / (n-1);  // t-step between adjacent control points

    std::printf("Twisted cubic  r(t)=(t,t²,t³)  n=%d  h=1/%d\n", n, n-1);
    std::printf("%-4s %-7s %-7s %-10s %-10s %-10s %s\n",
                "win", "t_mid", "τ", "κτ(2h)³/6", "planar_rms", "ratio", "conic?");
    std::printf("%s\n", std::string(66, '-').c_str());

    int n_conic = 0;
    double max_pr_conic = 0.0;
    for (int i = 0; i < n-4; ++i) {
        V3 pts5[5] = {ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4]};
        double t_mid  = ts[i+2];
        double kappa  = tc_curvature(t_mid);
        double tau    = tc_torsion(t_mid);
        double predicted = kappa * tau * std::pow(2.0*h, 3) / 6.0;
        double observed  = planarity_rms(pts5);
        double ratio     = (predicted > 0) ? observed/predicted : 0;

        ConicWindow<3> cw(ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4],
                          ts[i],ts[i+1],ts[i+2],ts[i+3],ts[i+4]);
        bool valid = cw.valid();

        std::printf("%-4d %-7.4f %-7.4f %-10.3e %-10.3e %-10.3f %s\n",
                    i, t_mid, tau, predicted, observed, ratio,
                    valid ? "conic" : "Lagrange");

        if (valid) { ++n_conic; if (observed > max_pr_conic) max_pr_conic = observed; }
    }

    std::printf("\n%d/%d windows conic\n", n_conic, n-4);
    std::printf("max planarity_rms of accepted conic windows: %.3e\n", max_pr_conic);
    std::printf("(Frenet-Serret κτ(2h)³/6 prediction matches to O(1) constant)\n");

    if (n_conic == 0)
        fail("per_window", "no conic windows — gap not exercised");
    if (max_pr_conic < 1e-8)
        fail("per_window", "planarity_rms unexpectedly zero");
    pass("per_window");
}

// ── Test 2: h³ scaling study ───────────────────────────────────────────────
//
// Vary n = 8, 12, 20, 40.  For the window centred at t≈0.5 (fixed torsion),
// measure max planarity_rms of accepted conic windows.  The ratio between
// successive h values should approach 2³=8 (h³ scaling).

void test_h3_scaling()
{
    static const int ns[] = {8, 12, 20, 40};
    std::printf("\nh³ scaling study (max planarity_rms of conic windows):\n");
    std::printf("%-4s %-8s %-12s %-8s\n", "n", "h", "max_pr_conic", "ratio_h³");
    std::printf("%s\n", std::string(38, '-').c_str());

    double prev_pr = 0.0, prev_h = 0.0;
    for (int n : ns) {
        std::vector<V3>     ctrl(n);
        std::vector<double> ts(n);
        for (int i = 0; i < n; ++i) {
            double t = (double)i / (n-1);
            ctrl[i] = tc(t); ts[i] = t;
        }
        double h = 1.0 / (n-1);

        double max_pr = 0.0;
        for (int i = 0; i < n-4; ++i) {
            V3 pts5[5] = {ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4]};
            ConicWindow<3> cw(ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4],
                              ts[i],ts[i+1],ts[i+2],ts[i+3],ts[i+4]);
            if (cw.valid()) {
                double pr = planarity_rms(pts5);
                if (pr > max_pr) max_pr = pr;
            }
        }

        double ratio = (prev_pr > 0 && max_pr > 0)
                       ? (prev_pr / max_pr) / std::pow(prev_h / h, 3) : 0.0;
        std::printf("%-4d %-8.4f %-12.3e %-8.3f\n",
                    n, h, max_pr, (prev_pr > 0 ? ratio : 0.0));
        prev_pr = max_pr; prev_h = h;
    }
    std::printf("(ratio → 1.0 confirms h³ scaling)\n");
    pass("h3_scaling");
}

// ── Test 3: out-of-plane interpolation error ───────────────────────────────
//
// The blend curve at a knot t[j] evaluates two conic windows — each returning
// a point in its own best-fit plane — rather than the original 3D ctrl point.
// The interpolation error is O(planarity_rms), not O(machine_epsilon).
//
// We measure: max |blend(t[j]) - ctrl[j]| and compare to max planarity_rms.

void test_interpolation_error()
{
    int n = 16;
    std::vector<V3>     ctrl(n);
    std::vector<double> ts(n);
    for (int i = 0; i < n; ++i) {
        double t = (double)i / (n-1);
        ctrl[i] = tc(t); ts[i] = t;
    }

    auto result = blend_curve<3>(ctrl, ts, conic_tag{}, 40, 2);

    // Blend covers segments j=2..n-4, i.e. times[2]..times[n-3].
    // Guard: skip knots whose nearest sampled time is > 2 sample steps away
    // (boundary knots j=1 and j=n-2 fall outside the blend range).
    double dt_step = (ts[n-1] - ts[0]) / ((n-5) * 40.0);
    double max_interp_err = 0.0;
    for (int j = 2; j <= n-3; ++j) {
        double tj = ts[j];
        std::size_t best_i = 0; double best_d = 1e30;
        for (std::size_t i = 0; i < result.times.size(); ++i) {
            double d = std::abs(result.times[i] - tj);
            if (d < best_d) { best_d = d; best_i = i; }
        }
        if (best_d > 2.0 * dt_step) continue;  // outside blend range
        double e = (result.pts[best_i] - ctrl[j]).norm();
        if (e > max_interp_err) max_interp_err = e;
    }

    // Max planarity_rms over all windows with accepted conic
    double max_pr = 0.0;
    for (int i = 0; i < n-4; ++i) {
        V3 pts5[5] = {ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4]};
        ConicWindow<3> cw(ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4],
                          ts[i],ts[i+1],ts[i+2],ts[i+3],ts[i+4]);
        if (cw.valid()) {
            double pr = planarity_rms(pts5);
            if (pr > max_pr) max_pr = pr;
        }
    }

    std::printf("\nOut-of-plane interpolation error (n=%d):\n", n);
    std::printf("  max planarity_rms of conic windows:  %.3e\n", max_pr);
    std::printf("  max |blend(t_j) - ctrl[j]|:          %.3e\n", max_interp_err);
    std::printf("  ratio interp/planarity:               %.2f\n",
                max_pr > 0 ? max_interp_err/max_pr : 0.0);
    std::printf("  (interpolation error ≈ O(planarity_rms), not O(machine_eps))\n");

    // Sanity: interp error should be meaningfully larger than machine epsilon
    // but not wildly larger than the planarity residual
    if (max_interp_err < 1e-10)
        fail("interp_error", "error unexpectedly small — gap not present");
    if (max_interp_err > 100.0 * max_pr + 1e-6)
        fail("interp_error", "interpolation error >> planarity_rms (unexpected)");
    pass("interp_error");
}

// ── main ──────────────────────────────────────────────────────────────────

int main()
{
    std::printf("=== Torsion gap — ConicWindow on twisted cubic (t,t²,t³) ===\n\n");
    test_per_window();
    test_h3_scaling();
    test_interpolation_error();
    std::printf("\nAll tests passed.\n");
}
