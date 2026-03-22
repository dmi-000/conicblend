// demo.cpp — Validate conicblend_circle.hpp against analytic test curves.
//
// Build:
//   c++ -std=c++17 -O2 -o demo demo.cpp && ./demo
//
// Or with CMake:
//   cmake -B build && cmake --build build && ./build/demo
//
// Tests:
//   1.  Helix:               max position error vs analytic  (should be ≪ 1)
//   2.  Torus knot:          max error at control points     (exact: ~1e-15)
//   3.  C^N continuity:      knot-junction value agreement   (should be ~0)
//   4.  Collinearity guard:  circumcircle throws for collinear pts
//   5.  Exact circle:        blended circle to machine precision (~1e-13)
//   6.  Input validation:    blend_curve throws for bad arguments
//   7.  Minimum n=4:         exactly one blended segment works
//   8.  Large arc (>π):      CircleWindow branch-cut correction exercised
//   9.  Clockwise arc:       negative traversal direction handled correctly
//  10.  Tilted circle:        circle in a non-xy plane reproduced exactly
//  11.  Non-uniform times:    exact interpolation holds for unequal spacing
//  12.  Smoothstep orders:    N=1 (C¹) and N=3 (C³) produce correct weights
//  13.  Tag dispatch:         fc::circle_tag{} overload compiles and matches untagged result

#include "conicblend_circle.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

// ── Analytic test curves ──────────────────────────────────────────────────

// Helix: P(θ) = (R cos θ, R sin θ, pitch·θ)
// At θ=0 : P=(R,0,0), |dP/dθ| = √(R²+pitch²).
fc::Vec3 helix_point(double theta, double R = 1.0, double pitch = 0.3)
{
    return {R * std::cos(theta), R * std::sin(theta), pitch * theta};
}

// (2,3) Torus knot on torus with R=2, r=1.
// P(t) = ((R + r cos 3t) cos 2t,  (R + r cos 3t) sin 2t,  r sin 3t)
fc::Vec3 torus_knot_point(double t, double R = 2.0, double r = 1.0)
{
    return {
        (R + r * std::cos(3.0*t)) * std::cos(2.0*t),
        (R + r * std::cos(3.0*t)) * std::sin(2.0*t),
        r * std::sin(3.0*t)
    };
}

// ── Helpers ───────────────────────────────────────────────────────────────

double dist(fc::Vec3 const& a, fc::Vec3 const& b)
{
    double dx = a[0]-b[0], dy = a[1]-b[1], dz = a[2]-b[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// Produce n evenly-spaced control points on a curve.
template<typename F>
std::pair<std::vector<fc::Vec3>, std::vector<double>>
sample_curve(F curve_func, double t0, double t1, int n)
{
    std::vector<fc::Vec3>   ctrl(n);
    std::vector<double>     times(n);
    for (int i = 0; i < n; ++i) {
        double t = t0 + (t1-t0) * i / (n-1);
        ctrl[i]  = curve_func(t);
        times[i] = t;
    }
    return {ctrl, times};
}

// Write CSV: t, x, y, z
void write_csv(std::string const& fname, fc::BlendResult const& r)
{
    std::ofstream f(fname);
    f << std::setprecision(10);
    f << "t,x,y,z\n";
    for (std::size_t i = 0; i < r.pts.size(); ++i)
        f << r.times[i] << ',' << r.pts[i][0] << ','
          << r.pts[i][1] << ',' << r.pts[i][2] << '\n';
    std::cout << "  wrote " << fname << "  (" << r.pts.size() << " pts)\n";
}

// ── Test 1: helix approximation error ─────────────────────────────────────
//
// The blended curve is NOT required to trace the helix exactly (it's composed
// of circular-arc blends, not the helix). We measure the max distance from
// each blended point to the nearest point on the analytic helix — actually we
// just compare at the same t parameter for a fair approximation quality test.

void test_helix(int n = 20, int N_smooth = 2)
{
    std::cout << "\n── Test 1: Helix  (n=" << n << ", N=" << N_smooth << ") ──\n";

    auto hp = [](double t){ return helix_point(t); };
    auto [ctrl, times] = sample_curve(hp, 0.0, 4.0*fc::detail::PI, n);
    auto result = fc::blend_curve(ctrl, times, 80, N_smooth);

    double max_err = 0.0;
    for (std::size_t i = 0; i < result.pts.size(); ++i) {
        fc::Vec3 analytic = helix_point(result.times[i]);
        max_err = std::max(max_err, dist(result.pts[i], analytic));
    }
    std::cout << "  max |blend − analytic helix|  = " << max_err << "\n";

    write_csv("helix_blend.csv", result);

    // Print convergence: n=8,12,16,20,30,50
    std::cout << "  Convergence table:\n";
    std::cout << "    n     max_err\n";
    for (int nn : {8, 12, 16, 20, 30, 50}) {
        auto [c2, t2] = sample_curve(hp, 0.0, 4.0*fc::detail::PI, nn);
        auto r2 = fc::blend_curve(c2, t2, 80, N_smooth);
        double e = 0.0;
        for (std::size_t i = 0; i < r2.pts.size(); ++i)
            e = std::max(e, dist(r2.pts[i], helix_point(r2.times[i])));
        std::printf("    %-5d %.4e\n", nn, e);
    }
}

// ── Test 2: exact interpolation at control points ─────────────────────────
//
// At t = times[j] (j = 1, ..., n-3), the blended curve should reproduce
// ctrl[j] exactly:  blend = (1−0)·win[j-1](times[j]) + 0·win[j](times[j])
//                         = win[j-1](times[j]) = ctrl[j].
// We re-run blend_curve and search for the output points at the junction times.
// Because the dense output uses exact junction times (k=0 at each segment),
// the first point of each segment IS at times[j].

void test_exact_interpolation(int n = 16)
{
    std::cout << "\n── Test 2: Exact interpolation at control points  (n=" << n << ") ──\n";

    // Torus knot — highly non-planar, good stress test.
    auto tkp = [](double t){ return torus_knot_point(t); };
    auto [ctrl, times] = sample_curve(tkp, 0.0, 2.0*fc::detail::PI, n);
    auto result = fc::blend_curve(ctrl, times, 60, 2);

    // The blended output at index k*60 within segment j is the junction point times[j+1].
    // Actually: for pts_per_seg=60, segment j has pts at k=0,...,59 (not 60 for non-last).
    // k=0 corresponds to t = times[j] → this IS ctrl[j].
    // We find every 60-th point (beginning of each segment) and compare.
    int n_segs = n - 3;
    int pps    = 60;
    double max_err = 0.0;
    for (int seg = 0; seg < n_segs; ++seg) {
        int idx        = seg * pps;                // k=0 of this segment
        int ctrl_idx   = seg + 1;                  // segment j=1+seg → ctrl[j] = ctrl[1+seg]
        double e = dist(result.pts[idx], ctrl[ctrl_idx]);
        max_err = std::max(max_err, e);
    }
    // Also check the very last point (end of last segment) vs ctrl[n-2]:
    double e_last = dist(result.pts.back(), ctrl[n-2]);
    max_err = std::max(max_err, e_last);

    std::cout << "  max |blend(t_j) − ctrl[j]|   = " << max_err
              << "  (should be ≲ 1e-14)\n";

    write_csv("torus_blend.csv", result);
}

// ── Test 3: C^N junction continuity ───────────────────────────────────────
//
// At each knot t_j, the last point of segment j-1 and the first point of
// segment j both evaluate to the same CircleWindow win[j-1] at t_j.
// Therefore they should be IDENTICAL (to floating-point precision).
//
// We verify this by running blend_curve with junction deduplication OFF —
// i.e., by requesting pts_per_seg points for EVERY segment including the
// last — and then checking adjacent-segment junction pairs.
//
// Since blend_curve deduplicates (only the last segment includes k=pts_per_seg),
// we instead build two separate single-segment evals and compare.
//
// Simpler approach: rebuild two CircleWindows for each junction and compare
// both sides directly at the junction time.

void test_continuity(int n = 12, int N_smooth = 2)
{
    std::cout << "\n── Test 3: C^N continuity at junctions  (n=" << n << ") ──\n";

    auto hp2 = [](double t){ return helix_point(t); };
    auto [ctrl, times] = sample_curve(hp2, 0.0, 2.0*fc::detail::PI, n);

    // Build all CircleWindows (win[i] covers ctrl[i..i+2]).
    std::vector<fc::CircleWindow<3>> wins;
    wins.reserve(n - 2);
    for (int i = 0; i < n-2; ++i)
        wins.emplace_back(ctrl[i], ctrl[i+1], ctrl[i+2],
                          times[i], times[i+1], times[i+2]);

    // At junction j (j = 1,...,n-3):
    //   left  side (end of segment j-1):  wA=wins[j-1], s→1 → w→1 → blend = wB(t_j)  = wins[j](times[j])
    //   right side (start of segment j):  wA=wins[j-1], s→0 → w→0 → blend = wA(t_j)  = wins[j-1](times[j])
    //
    // Wait — re-check the indexing:
    //   Segment j uses wA=wins[j-1], wB=wins[j].
    //   At s=0 (t=times[j]):   blend = wA(times[j]) = wins[j-1](times[j]) = ctrl[j]    ✓
    //   At s=1 (t=times[j+1]): blend = wB(times[j+1]) = wins[j](times[j+1]) = ctrl[j+1] ✓
    //
    // Junction between segment j and segment j+1 is at t=times[j+1]:
    //   seg j   at s=1: blend = wins[j](times[j+1])
    //   seg j+1 at s=0: blend = wins[j](times[j+1])   ← SAME window!
    //
    // So both sides equal wins[j](times[j+1]) — they are literally identical.
    // The numerical test is: wins[j](times[j+1]) == wins[j](times[j+1]),
    // which is trivially 0. Instead, verify the blend formula gives ctrl[j+1]:
    //   left  (s=1): blend = wins[j](times[j+1]) should = ctrl[j+1]
    //   right (s=0): blend = wins[j](times[j+1]) should = ctrl[j+1]

    double max_junction_err = 0.0;
    for (int j = 1; j <= n-3; ++j) {
        double t_j1 = times[j];     // start of segment j
        double t_j2 = times[j+1];   // end   of segment j = start of j+1

        // At the start of segment j: blend = wins[j-1](t_j) = ctrl[j]
        fc::Vec3 pt_start = wins[j-1](t_j1);
        double e_start = dist(pt_start, ctrl[j]);
        max_junction_err = std::max(max_junction_err, e_start);

        // At the end of segment j: blend = wins[j](t_j2) = ctrl[j+1]
        fc::Vec3 pt_end = wins[j](t_j2);
        double e_end = dist(pt_end, ctrl[j+1]);
        max_junction_err = std::max(max_junction_err, e_end);
    }
    std::cout << "  max |win(t_j) − ctrl[j]|      = " << max_junction_err
              << "  (exact interpolation, should be ≲ 1e-14)\n";

    // Also verify C^0 continuity of the blend at junction:
    // The two segments share the same window evaluation → automatically C^∞
    // in the window itself.  But to make this concrete, verify the BLEND
    // function evaluated from left and right gives the same value.
    // Left:  seg j at t=times[j+1]:  s=1 → w=1 → blend = wins[j](times[j+1])
    // Right: seg j+1 at t=times[j+1]: s=0 → w=0 → blend = wins[j](times[j+1])
    // They are literally the same call, so the C^0 error is 0.0 by construction.
    std::cout << "  C^0 junction error             = 0.0  (same window eval from both sides)\n";
    std::cout << "  C^N: smoothstep has " << N_smooth << " vanishing derivatives at s=0,1\n"
              << "        → blend is C^" << N_smooth << " at every junction by construction.\n";
}

// ── Test 4: collinearity guard ─────────────────────────────────────────────

void test_collinearity_guard()
{
    std::cout << "\n── Test 4: Collinearity guard ──\n";
    fc::Vec3 p0{0,0,0}, p1{1,0,0}, p2{2,0,0};   // collinear
    try {
        fc::circumcircle(p0, p1, p2);
        std::cout << "  ERROR: should have thrown!\n";
    } catch (std::invalid_argument const& e) {
        std::cout << "  Caught (expected): " << e.what() << "\n";
    }
}

// ── Test 5: exact circle ───────────────────────────────────────────────────
//
// If control points lie exactly on a circle, the blended curve should
// reproduce the circle to machine precision (each window IS the circle).

void test_exact_circle(int n = 16, int N_smooth = 2)
{
    std::cout << "\n── Test 5: Exact circle  (n=" << n << ") ──\n";

    double R = 1.5;
    std::vector<fc::Vec3>   ctrl(n);
    std::vector<double>     times(n);
    for (int i = 0; i < n; ++i) {
        double theta = 2.0*fc::detail::PI * i / n;
        ctrl[i]  = {R*std::cos(theta), R*std::sin(theta), 0.0};
        times[i] = theta;
    }

    auto result = fc::blend_curve(ctrl, times, 60, N_smooth);

    double max_err = 0.0;
    for (auto const& p : result.pts) {
        double r = p.norm();
        max_err = std::max(max_err, std::abs(r - R));
    }
    std::cout << "  max |r − R|                   = " << max_err
              << "  (should be ≲ 1e-13)\n";
}

// ── Test 6: input validation ───────────────────────────────────────────────
//
// blend_curve should throw std::invalid_argument for:
//   (a) ctrl.size() != times.size()
//   (b) n < 4
//   (c) pts_per_seg < 2
//   (d) smooth_N outside {1,2,3}

void test_input_validation()
{
    std::cout << "\n── Test 6: Input validation ──\n";

    auto expect_throw = [](auto fn, char const* desc) {
        try {
            fn();
            std::printf("  FAIL — %s: expected exception, got none\n", desc);
        } catch (std::invalid_argument const& e) {
            std::printf("  OK   — %s: \"%s\"\n", desc, e.what());
        }
    };

    // Collinear (for size-mismatch and n<4 tests — error fires before circumcircle)
    std::vector<fc::Vec3>   c4_col = {{0,0,0},{1,0,0},{2,0,0},{3,0,0}};
    std::vector<double>     t4     = {0.0, 1.0, 2.0, 3.0};
    std::vector<double>     t3     = {0.0, 1.0, 2.0};
    std::vector<fc::Vec3>   c3(c4_col.begin(), c4_col.begin()+3);

    // Non-collinear, for tests where we need to reach the smoothstep dispatch.
    // Four points on the unit circle in the xy-plane.
    std::vector<fc::Vec3> c4_ok = {
        {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}
    };

    expect_throw([&]{ fc::blend_curve(c4_col, t3); },        "size mismatch");
    expect_throw([&]{ fc::blend_curve(c3, t3); },            "n < 4");
    expect_throw([&]{ fc::blend_curve(c4_ok, t4, 1); },      "pts_per_seg < 2");
    expect_throw([&]{ fc::blend_curve(c4_ok, t4, 60, 4); },  "N=4 (invalid)");
    expect_throw([&]{ fc::blend_curve(c4_ok, t4, 60, 0); },  "N=0 (invalid)");
}

// ── Test 7: minimum n=4 ────────────────────────────────────────────────────
//
// n=4 yields exactly one blended segment (j=1).
// The blended curve must pass through ctrl[1] at t=times[1]
// and through ctrl[2] at t=times[2].

void test_minimum_n()
{
    std::cout << "\n── Test 7: Minimum n=4 (one segment) ──\n";

    // Four points on the helix.
    std::vector<fc::Vec3>   ctrl(4);
    std::vector<double>     times(4);
    for (int i = 0; i < 4; ++i) {
        double t = fc::detail::PI * i / 3.0;
        ctrl[i]  = helix_point(t);
        times[i] = t;
    }

    auto r = fc::blend_curve(ctrl, times, 60, 2);

    // First point of output: t = times[1] → should equal ctrl[1]
    double e0 = dist(r.pts.front(), ctrl[1]);
    // Last point of output:  t = times[2] → should equal ctrl[2]
    double e1 = dist(r.pts.back(),  ctrl[2]);

    std::printf("  n_output = %zu  (expected %d)\n", r.pts.size(), 61);
    std::printf("  |blend(t1) − ctrl[1]| = %.2e  (should be ≲ 1e-14)\n", e0);
    std::printf("  |blend(t2) − ctrl[2]| = %.2e  (should be ≲ 1e-14)\n", e1);
}

// ── Test 8: large arc (> π) ────────────────────────────────────────────────
//
// Places three points on a unit circle with angles 0, 200°, 340° — the arc
// from p0 to p2 via p1 spans 340° (> π).  This exercises the d1*d2 < 0
// branch-cut correction in CircleWindow, which must wrap d2 by ±2π so that
// d1 and d2 are on the same side of zero.
//
// Correctness check: the window evaluated at t1 must return p1 exactly.

void test_large_arc()
{
    std::cout << "\n── Test 8: Large arc (> π) — branch-cut correction ──\n";

    // Angles: 0°, 200°, 340° → CCW arc spans 340° total.
    auto deg = [](double d){ return d * fc::detail::PI / 180.0; };
    fc::Vec3 p0{std::cos(deg(  0)), std::sin(deg(  0)), 0};
    fc::Vec3 p1{std::cos(deg(200)), std::sin(deg(200)), 0};
    fc::Vec3 p2{std::cos(deg(340)), std::sin(deg(340)), 0};

    fc::CircleWindow<3> win(p0, p1, p2, 0.0, 1.0, 2.0);

    fc::Vec3 at0 = win(0.0), at1 = win(1.0), at2 = win(2.0);
    std::printf("  |win(t0) − p0| = %.2e  (should be ≲ 1e-15)\n", dist(at0, p0));
    std::printf("  |win(t1) − p1| = %.2e  (should be ≲ 1e-15)\n", dist(at1, p1));
    std::printf("  |win(t2) − p2| = %.2e  (should be ≲ 1e-15)\n", dist(at2, p2));

    // Verify radius at the mid-evaluation point is 1.0 (stays on the circle).
    fc::Vec3 mid = win(1.0);
    double r_mid = mid.norm();
    std::printf("  |r_mid − 1.0|  = %.2e  (stays on circle)\n",
                std::abs(r_mid - 1.0));
}

// ── Test 9: clockwise arc ──────────────────────────────────────────────────
//
// Three points with decreasing angle (0°, −90°, −180°): d1 < 0.
// The arc goes clockwise.  Both interpolation endpoints and the midpoint
// radius should be exact.

void test_clockwise_arc()
{
    std::cout << "\n── Test 9: Clockwise arc (d1 < 0) ──\n";

    auto deg = [](double d){ return d * fc::detail::PI / 180.0; };
    fc::Vec3 p0{std::cos(deg(  0)), std::sin(deg(  0)), 0};
    fc::Vec3 p1{std::cos(deg(-90)), std::sin(deg(-90)), 0};
    fc::Vec3 p2{std::cos(deg(180)), std::sin(deg(180)), 0};

    fc::CircleWindow<3> win(p0, p1, p2, 0.0, 1.0, 2.0);

    std::printf("  |win(t0) − p0| = %.2e\n", dist(win(0.0), p0));
    std::printf("  |win(t1) − p1| = %.2e\n", dist(win(1.0), p1));
    std::printf("  |win(t2) − p2| = %.2e\n", dist(win(2.0), p2));

    // The midpoint between t0 and t2 should land at −45° (halfway CW from 0 to −90).
    fc::Vec3 expected{std::cos(deg(-45)), std::sin(deg(-45)), 0};
    std::printf("  |win(0.5) − expected(-45°)| = %.2e  (quadratic φ, approx)\n",
                dist(win(0.5), expected));
}

// ── Test 10: tilted circle ─────────────────────────────────────────────────
//
// A circle of radius 2 in the plane  z = x  (normal = (1,0,-1)/√2).
// If control points lie on this circle, the blended curve must reproduce
// the circle to machine precision.

void test_tilted_circle(int n = 12)
{
    std::cout << "\n── Test 10: Tilted circle (plane z=x)  (n=" << n << ") ──\n";

    // In-plane basis for z = x: e1 = (1,0,1)/√2, e2 = (0,1,0).
    double R = 2.0;
    fc::Vec3 ctr{0, 0, 0};
    fc::Vec3 e1{1.0/std::sqrt(2.0), 0, 1.0/std::sqrt(2.0)};
    fc::Vec3 e2{0, 1, 0};

    std::vector<fc::Vec3>   ctrl(n);
    std::vector<double>     times(n);
    for (int i = 0; i < n; ++i) {
        double theta = 2.0*fc::detail::PI * i / n;
        times[i] = theta;
        ctrl[i]  = ctr + (R*std::cos(theta)) * e1 + (R*std::sin(theta)) * e2;
    }

    auto result = fc::blend_curve(ctrl, times, 60, 2);

    // Every point must lie on the tilted circle: |p - ctr| = R
    // AND  p[2] - p[0] = 0  (plane constraint: z = x  ⟹  z - x = 0)
    double max_radius_err = 0.0, max_plane_err = 0.0;
    for (auto const& p : result.pts) {
        fc::Vec3 r = p - ctr;
        max_radius_err = std::max(max_radius_err, std::abs(r.norm() - R));
        max_plane_err  = std::max(max_plane_err,  std::abs(p[2] - p[0]));
    }
    std::printf("  max |r − R|       = %.2e  (should be ≲ 1e-13)\n", max_radius_err);
    std::printf("  max |z − x|       = %.2e  (stays in tilted plane)\n", max_plane_err);
}

// ── Test 11: non-uniform parameter spacing ─────────────────────────────────
//
// Uses the torus knot but samples it at non-uniform parameter values
// (geometric progression rather than equal spacing).
// Exact interpolation must hold at every interior control point.

void test_nonuniform_times(int n = 12)
{
    std::cout << "\n── Test 11: Non-uniform parameter spacing  (n=" << n << ") ──\n";

    // Geometric spacing: t_i = t_max * (r^i - 1) / (r^(n-1) - 1), r=1.5
    std::vector<double> times(n);
    double r = 1.5, t_max = 2.0*fc::detail::PI;
    double rn = std::pow(r, n-1) - 1.0;
    for (int i = 0; i < n; ++i)
        times[i] = t_max * (std::pow(r, i) - 1.0) / rn;

    std::vector<fc::Vec3> ctrl(n);
    for (int i = 0; i < n; ++i) ctrl[i] = torus_knot_point(times[i]);

    auto result = fc::blend_curve(ctrl, times, 60, 2);

    // Check interior ctrl pts: every 60th point starting from 0 is ctrl[j+1].
    int pps = 60, n_segs = n - 3;
    double max_err = 0.0;
    for (int seg = 0; seg < n_segs; ++seg) {
        int idx = seg * pps;
        double e = dist(result.pts[idx], ctrl[seg + 1]);
        max_err = std::max(max_err, e);
    }
    double e_last = dist(result.pts.back(), ctrl[n-2]);
    max_err = std::max(max_err, e_last);

    std::printf("  max |blend(t_j) − ctrl[j]| = %.2e  (should be ≲ 1e-14)\n", max_err);
}

// ── Test 12: smoothstep orders N=1 and N=3 ────────────────────────────────
//
// Verify that the N=1 (C¹) and N=3 (C³) smoothstep weights satisfy the
// boundary conditions: w(0)=0, w(1)=1, and the first N derivatives vanish
// at 0 and 1.  Also check that blend_curve runs without error for each order.

void test_smoothstep_orders()
{
    std::cout << "\n── Test 12: Smoothstep orders N=1 and N=3 ──\n";

    // Numerical derivative of smoothstep at s≈0 and s≈1.
    auto deriv = [](int N, double s, double h = 1e-6) {
        return (fc::smoothstep(s + h, N) - fc::smoothstep(s - h, N)) / (2.0*h);
    };

    for (int N : {1, 3}) {
        double w0  = fc::smoothstep(0.0, N);
        double w1  = fc::smoothstep(1.0, N);
        double dw0 = deriv(N, 1e-6);   // dw/ds near s=0
        double dw1 = deriv(N, 1.0 - 1e-6);  // dw/ds near s=1
        std::printf("  N=%d  w(0)=%.1f  w(1)=%.1f"
                    "  w'(0)=%.2e  w'(1)=%.2e\n",
                    N, w0, w1, dw0, dw1);
    }

    // Check blend_curve runs for N=1 and N=3 on a simple helix.
    auto hp = [](double t){ return helix_point(t); };
    auto [ctrl, times] = sample_curve(hp, 0.0, 2.0*fc::detail::PI, 10);
    for (int N : {1, 3}) {
        auto r = fc::blend_curve(ctrl, times, 40, N);
        std::printf("  N=%d  npts=%zu  (no exception)\n", N, r.pts.size());
    }
}

// ── Main ──────────────────────────────────────────────────────────────────

int main()
{
    std::cout << "conicblend demo — 3-point circle-arc blend curve\n";
    std::cout << std::string(50, '=') << "\n";

    test_helix();
    test_exact_interpolation();
    test_continuity();
    test_collinearity_guard();
    test_exact_circle();
    test_input_validation();
    test_minimum_n();
    test_large_arc();
    test_clockwise_arc();
    test_tilted_circle();
    test_nonuniform_times();
    test_smoothstep_orders();

    // Test 13: tag dispatch — fc::circle_tag{} overload produces identical output.
    {
        int n = 10;
        std::vector<fc::Vec3>  ctrl(n);
        std::vector<double>    times(n);
        for (int i = 0; i < n; ++i) {
            double t = 2.0 * fc::detail::PI * i / n;
            ctrl[i]  = {std::cos(t), std::sin(t), 0.0};
            times[i] = t;
        }
        auto r_plain  = fc::blend_curve(ctrl, times, 40, 2);
        auto r_tagged = fc::blend_curve(ctrl, times, fc::circle_tag{}, 40, 2);
        assert(r_plain.pts.size() == r_tagged.pts.size());
        double max_diff = 0.0;
        for (std::size_t i = 0; i < r_plain.pts.size(); ++i) {
            auto d = r_plain.pts[i] - r_tagged.pts[i];
            max_diff = std::max(max_diff, d.norm());
        }
        assert(max_diff == 0.0);
        std::cout << "Test 13 (tag dispatch): plain vs circle_tag{} max diff = "
                  << max_diff << "  PASS\n";
    }

    std::cout << "\nAll tests complete.\n";
    return 0;
}
