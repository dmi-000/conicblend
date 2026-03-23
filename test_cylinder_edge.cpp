// test_cylinder_edge.cpp — edge-case tests for ConicWindow + CylinderWindow
//
// Tests new code paths added 2026-03-23:
//   line_mode_ + fit_error (ConicWindow), cyl_best_fit + used_best_fit (CylinderWindow)
//
// Compile:
//   g++ -std=c++17 -O2 -I. \
//       -isystem /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 \
//       -o test_cylinder_edge test_cylinder_edge.cpp

#include "conicblend_cylinder.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstring>

using namespace fc;
using V2 = VecN<2>;
using V3 = VecN<3>;

static int n_pass = 0, n_fail = 0;

#define CHECK(desc, cond) do { \
    bool _ok = (cond); \
    std::printf("%s  %s\n", _ok ? "PASS" : "FAIL", (desc)); \
    if (_ok) ++n_pass; else ++n_fail; \
} while(0)

// ── helpers ──────────────────────────────────────────────────────────────────

static V3 helix5_pts[5];
static double helix5_ts[5];

static void make_helix5() {
    for (int k = 0; k < 5; ++k) {
        helix5_ts[k] = 0.1 * k;
        double t = helix5_ts[k];
        helix5_pts[k] = V3{std::cos(2*t), std::sin(2*t), 0.5*t};
    }
}

static double win_scale(const V3 pts[5]) {
    double ws = 0;
    for (int i = 0; i < 5; ++i)
        for (int j = i+1; j < 5; ++j)
            ws = std::max(ws, (pts[i]-pts[j]).norm());
    return ws;
}

// ── Test 1: ConicWindow<2> collinear → line mode ──────────────────────────
static void test1_line_mode() {
    std::printf("\n── T1: ConicWindow<2> collinear input → line mode ──────────────────\n");

    V2 pts[5] = {{-2,0},{-1,0},{0,0},{1,0},{2,0}};
    double ts[5] = {0,1,2,3,4};
    ConicWindow<2> w(pts[0],pts[1],pts[2],pts[3],pts[4], ts[0],ts[1],ts[2],ts[3],ts[4]);

    CHECK("T1a valid=true",          w.valid());
    CHECK("T1b fit_error=0",         w.fit_error() == 0.0);

    double max_err = 0;
    for (int k = 0; k < 5; ++k) {
        V2 p = w(ts[k]);
        max_err = std::max(max_err, (p - pts[k]).norm());
    }
    std::printf("     max knot error = %.2e\n", max_err);
    CHECK("T1c knots exact (< 1e-14)", max_err < 1e-14);
}

// ── Test 2: fit_error = 0 in line mode at large scale ─────────────────────
//
// 5 collinear points on the oblique line y = 2x + 1 at two different scales.
// Both windows use line mode (Lagrange, fit_error=0 exactly).
// Confirms the dimensionless win_scale normalisation does not overflow/underflow
// and does not spuriously reject large-coordinate inputs.
static void test2_scale_invariance() {
    std::printf("\n── T2: fit_error = 0 in line mode at unit and large scale ──────────\n");

    double ts[5] = {0,1,2,3,4};
    V2 pts[5], ptsK[5];
    for (int k = 0; k < 5; ++k) {
        double x = (double)k + 1.0;
        pts[k]  = V2{x,       2.0*x + 1.0};         // scale 1
        ptsK[k] = V2{x*1e6,  (2.0*x + 1.0)*1e6};   // scale 1e6
    }

    ConicWindow<2> w1(pts[0], pts[1], pts[2], pts[3], pts[4], ts[0],ts[1],ts[2],ts[3],ts[4]);
    ConicWindow<2> w2(ptsK[0],ptsK[1],ptsK[2],ptsK[3],ptsK[4],ts[0],ts[1],ts[2],ts[3],ts[4]);

    CHECK("T2a both valid (line mode, any scale)",  w1.valid() && w2.valid());
    std::printf("     fit_error unit=%.2e  large=%.2e\n",
                w1.fit_error(), w2.fit_error());
    CHECK("T2b fit_error = 0 at unit scale",         w1.fit_error() == 0.0);
    CHECK("T2c fit_error = 0 at large scale (1e6)",  w2.fit_error() == 0.0);
}

// ── Test 3: CylinderWindow circle (zero-pitch geodesic) ───────────────────
static void test3_circle() {
    std::printf("\n── T3: CylinderWindow circle (zero-pitch geodesic) ─────────────────\n");

    V3 pts[5];
    double ts[5];
    for (int k = 0; k < 5; ++k) {
        double a = M_PI/3.0 * k;  // 0°, 60°, 120°, 180°, 240°
        pts[k] = V3{std::cos(a), std::sin(a), 0.0};
        ts[k]  = (double)k;
    }
    CylinderWindow<3> w(pts[0],pts[1],pts[2],pts[3],pts[4], ts[0],ts[1],ts[2],ts[3],ts[4]);

    CHECK("T3a valid=true",          w.valid());
    CHECK("T3b no best-fit used",    !w.used_best_fit());

    double max_err = 0;
    for (int k = 0; k < 5; ++k) {
        V3 p = w(ts[k]);
        max_err = std::max(max_err, (p - pts[k]).norm());
    }
    std::printf("     max knot error = %.2e\n", max_err);
    CHECK("T3c knots exact (< 1e-12)", max_err < 1e-12);
}

// ── Test 4: CylinderWindow helix knot interpolation ───────────────────────
static void test4_helix_knots() {
    std::printf("\n── T4: CylinderWindow helix knot interpolation ──────────────────────\n");

    make_helix5();
    CylinderWindow<3> w(
        helix5_pts[0],helix5_pts[1],helix5_pts[2],helix5_pts[3],helix5_pts[4],
        helix5_ts[0], helix5_ts[1], helix5_ts[2], helix5_ts[3], helix5_ts[4]);

    CHECK("T4a valid=true",          w.valid());
    CHECK("T4b no best-fit used",    !w.used_best_fit());

    double max_err = 0;
    for (int k = 0; k < 5; ++k) {
        V3 p = w(helix5_ts[k]);
        max_err = std::max(max_err, (p - helix5_pts[k]).norm());
    }
    std::printf("     max knot error = %.2e\n", max_err);
    CHECK("T4c knots exact (< 1e-12)", max_err < 1e-12);
}

// ── Test 5: CylinderWindow on known cylinder axis=(1,0,0) r=1 ─────────────
//
// pts[k] = (k, cos(k), sin(k)) for k=0..4.  All lie on the cylinder with
// axis=(1,0,0) through (0,0,0) and radius=1 (since y²+z²=cos²k+sin²k=1).
// cyl_solve finds (a,c)=(0,0) → axis ∥ x, b=d=0 exactly.
// Unrolled: (r·φ, z) = (k, k) — collinear → ConicWindow<2> line mode.
// Re-roll should reproduce pts[k] exactly.
static void test5_x_axis_cylinder() {
    std::printf("\n── T5: CylinderWindow axis=(1,0,0) r=1, knots exact ────────────────\n");

    V3 pts[5];
    double ts[5];
    for (int k = 0; k < 5; ++k) {
        pts[k] = V3{(double)k, std::cos((double)k), std::sin((double)k)};
        ts[k]  = (double)k;
    }
    CylinderWindow<3> w(pts[0],pts[1],pts[2],pts[3],pts[4], ts[0],ts[1],ts[2],ts[3],ts[4]);

    CHECK("T5a axis=(1,0,0) cylinder: valid=true", w.valid());
    CHECK("T5b axis=(1,0,0) cylinder: no best-fit", !w.used_best_fit());

    double max_err = 0;
    for (int k = 0; k < 5; ++k) {
        V3 p = w(ts[k]);
        max_err = std::max(max_err, (p - pts[k]).norm());
    }
    std::printf("     max knot error = %.2e\n", max_err);
    CHECK("T5c axis=(1,0,0) cylinder: knots exact (< 1e-12)", max_err < 1e-12);
}

// ── Test 6: cyl_best_fit unit — helix → accurate axis ─────────────────────
static void test6_best_fit_axis() {
    std::printf("\n── T6: cyl_best_fit unit — helix axis recovery ──────────────────────\n");

    make_helix5();
    double ws = win_scale(helix5_pts);
    auto bf = detail::cyl_best_fit(helix5_pts, ws);

    double axis_dot = std::abs(bf.u_hat.dot(V3{0,0,1}));
    std::printf("     axis·(0,0,1) = %.6f   max_surf_res = %.2e\n",
                axis_dot, bf.max_surf_res);
    CHECK("T6a helix axis ≈ (0,0,1)  (dot > 0.9999)", axis_dot > 0.9999);
    CHECK("T6b helix residuals ≈ 0   (< 1e-6)",        bf.max_surf_res < 1e-6);
}

// ── Test 7: cyl_best_fit gate — non-cylindrical point set ─────────────────
//
// 5 points with no cylindrical structure:
//   (1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1)
// For any axis direction the perpendicular projections span very different
// radii → Coope LS gives large surface residual → gate fires.
static void test7_best_fit_gate() {
    std::printf("\n── T7: cyl_best_fit gate — non-cylindrical point set ───────────────\n");

    // 5 points: origin + 3 axis endpoints + body diagonal.
    // No right circular cylinder can pass close to all five:
    //   the origin and (1,1,1) project to distances 0 vs √(2/3) for axis=(1,1,1)/√3,
    //   while (2,0,0),(0,2,0),(0,0,2) project to distance 2√(2/3) ≈ 1.63.
    // The best-fit cylinder must accommodate this large spread → max_surf_res > 0.05.
    V3 pts[5] = {{0,0,0},{2,0,0},{0,2,0},{0,0,2},{1,1,1}};

    double ws = win_scale(pts);
    auto bf = detail::cyl_best_fit(pts, ws);

    std::printf("     max_surf_res = %.3f  (gate threshold = 0.05)\n", bf.max_surf_res);
    CHECK("T7  non-cylindrical set: max_surf_res > 5e-2", bf.max_surf_res > 5e-2);
}

// ── Test 9: CylinderWindow coplanar input → valid()=false ─────────────────
//
// 5 coplanar points (z=0 plane, not collinear, not on any finite cylinder).
// cyl_solve finds no exact cylinder; cyl_best_fit gets max_surf_res > 5e-2
// because no finite cylinder can approximate a flat configuration well.
// CylinderWindow must be invalid — blend_curve falls back to LagrangeWindow<3>.
//
// NOTE: best-fit cylinder path (used_best_fit()=true) is NOT currently exercised
// by any test.  Triggering it requires 5 points where every exact-fit cylinder
// produces an invalid ConicWindow<2> (non-monotone φ or large fit_error), yet
// the grid-search best-fit finds a cylinder with max_surf_res < 5e-2.  Such a
// configuration is hard to construct analytically; coverage of that path is left
// to future work or a dedicated empirical regression.
static void test9_coplanar() {
    std::printf("\n── T9: CylinderWindow coplanar input → valid()=false ───────────────\n");

    // 5 points in the z=0 plane, in general position (not collinear)
    V3 pts[5] = {{1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{0.5,0.5,0}};
    double ts[5] = {0,1,2,3,4};
    CylinderWindow<3> w(pts[0],pts[1],pts[2],pts[3],pts[4], ts[0],ts[1],ts[2],ts[3],ts[4]);

    CHECK("T9a coplanar: valid()=false",       !w.valid());
    CHECK("T9b coplanar: used_best_fit=false",  !w.used_best_fit());
}

// ── Test 10: CylinderWindow tilted cylinder axis=(1,1,0)/√2 ───────────────
//
// Explicitly-constructed geodesic on a cylinder whose axis is tilted at 45°
// in the xy-plane (not aligned with any coordinate axis).  Tests that cyl_solve
// finds non-axis-aligned cylinders correctly and the re-roll is exact.
//
// Construction: u=(1,1,0)/√2, v=(1,-1,0)/√2, w=(0,0,1), offset=(0,0,0).
// pts[k] = u*k + v*cos(k) + w*sin(k)  (parameter k=0..4).
// These lie on the cylinder with axis u, radius=1, centered at the origin.
// The unrolled coordinates are (φ=k, z=k) — collinear → ConicWindow<2> line mode.
static void test10_tilted_cylinder() {
    std::printf("\n── T10: CylinderWindow tilted axis=(1,1,0)/√2, r=1, knots exact ───\n");

    const double s2 = std::sqrt(2.0);
    V3 u_ax = V3{1/s2, 1/s2, 0};
    V3 v_ax = V3{1/s2,-1/s2, 0};
    V3 w_ax = V3{0,   0,     1};

    V3 pts[5];
    double ts[5];
    for (int k = 0; k < 5; ++k) {
        double t = (double)k;
        ts[k]  = t;
        pts[k] = u_ax*t + v_ax*std::cos(t) + w_ax*std::sin(t);
    }
    CylinderWindow<3> w(pts[0],pts[1],pts[2],pts[3],pts[4], ts[0],ts[1],ts[2],ts[3],ts[4]);

    CHECK("T10a tilted cyl: valid=true",      w.valid());
    CHECK("T10b tilted cyl: no best-fit",     !w.used_best_fit());

    double max_err = 0;
    for (int k = 0; k < 5; ++k) {
        V3 p = w(ts[k]);
        max_err = std::max(max_err, (p - pts[k]).norm());
    }
    std::printf("     max knot error = %.2e\n", max_err);
    CHECK("T10c tilted cyl: knots exact (< 1e-12)", max_err < 1e-12);
}

// ── Test 8: full blend regression — helix n=8 ─────────────────────────────
static void test8_blend_regression() {
    std::printf("\n── T8: full blend regression — helix n=8 ───────────────────────────\n");

    auto helix = [](double t) { return V3{std::cos(2*t), std::sin(2*t), 0.5*t}; };

    int n = 8;
    std::vector<V3>     ctrl(n);
    std::vector<double> times(n);
    for (int i = 0; i < n; ++i) {
        times[i] = 2*M_PI * i / (n-1);
        ctrl[i]  = helix(times[i]);
    }

    auto r = blend_curve(ctrl, times, cylinder_tag{}, 200, 2);

    double max_d = 0;
    for (std::size_t i = 0; i < r.pts.size(); ++i) {
        double e = (r.pts[i] - helix(r.times[i])).norm();
        max_d = std::max(max_d, e);
    }
    std::printf("     max_dev = %.2e\n", max_d);
    CHECK("T8  helix n=8 max_dev < 1e-14", max_d < 1e-14);
}

// ── main ──────────────────────────────────────────────────────────────────
int main()
{
    std::printf("=== test_cylinder_edge: 10 edge-case tests ===\n");

    test1_line_mode();
    test2_scale_invariance();
    test3_circle();
    test4_helix_knots();
    test5_x_axis_cylinder();
    test6_best_fit_axis();
    test7_best_fit_gate();
    test9_coplanar();
    test10_tilted_cylinder();
    test8_blend_regression();

    std::printf("\n%d/%d passed\n", n_pass, n_pass + n_fail);
    return n_fail ? 1 : 0;
}
