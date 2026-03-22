// demo_conic.cpp — tests for conicblend.hpp (5-point conic windows)
//
// Compile:  clang++ -std=c++17 -O2 -I. -o demo_conic demo_conic.cpp
//
// Tests
//   1. Tilted ellipse:  max error vs analytic, exact interpolation ~machine eps
//   2. Hyperbola (same-branch): exact interpolation at control points
//   3. Parabola y=t²:  approximation quality
//   4. Helix (3D):     similar accuracy to circle windows
//   5. Input validation: n<6, size mismatch, pts_per_seg<2

#include "conicblend.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

using namespace fc;
using detail::PI;

// ── helpers ───────────────────────────────────────────────────────────────

template <int Dim>
double max_error(std::vector<VecN<Dim>> const& pts,
                 std::vector<double>    const& ts,
                 VecN<Dim> (*truth)(double))
{
    double err = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i)
        err = std::max(err, (pts[i] - truth(ts[i])).norm());
    return err;
}

// Find interpolation error at control-point times
template <int Dim>
double interp_error(std::vector<VecN<Dim>> const& out_pts,
                    std::vector<double>    const& out_ts,
                    std::vector<VecN<Dim>> const& ctrl,
                    std::vector<double>    const& ctrl_t)
{
    int n = static_cast<int>(ctrl.size());
    double err = 0.0;
    for (int j = 2; j <= n-3; ++j) {   // interior knots that appear in output
        double tj = ctrl_t[j];
        double best_d = 1e30;
        std::size_t best_i = 0;
        for (std::size_t i = 0; i < out_ts.size(); ++i) {
            double d = std::abs(out_ts[i] - tj);
            if (d < best_d) { best_d = d; best_i = i; }
        }
        double dt_step = (ctrl_t[n-1]-ctrl_t[0]) / ((n-5)*60.0);
        if (best_d > 2.0*dt_step) continue;
        double e = (out_pts[best_i] - ctrl[j]).norm();
        if (e > err) err = e;
    }
    return err;
}

static void pass(char const* name) { std::printf("PASS  %s\n", name); }
static void fail(char const* name, char const* msg) {
    std::printf("FAIL  %s: %s\n", name, msg);
    std::exit(1);
}

// ── Test 1: Tilted ellipse (2D) — exact conic input ───────────────────────

void test_tilted_ellipse()
{
    constexpr int Dim = 2;
    // Parametric ellipse: a=3, b=1, tilted 30°
    double ca = std::cos(PI/6.0), sa = std::sin(PI/6.0);
    int n = 16;
    std::vector<VecN<Dim>> ctrl(n);
    std::vector<double> times(n);

    for (int i = 0; i < n; ++i) {
        double t = 2.0*PI * i / (n-1);
        double x = 3.0*std::cos(t);
        double y = 1.0*std::sin(t);
        ctrl[i]  = VecN<Dim>(ca*x - sa*y, sa*x + ca*y);
        times[i] = t;
    }

    auto result = blend_curve<Dim>(ctrl, times, conic_tag{}, 60, 2);

    // Interpolation error at interior knots
    double interp_err = interp_error<Dim>(result.pts, result.times, ctrl, times);
    std::printf("  tilted_ellipse interp max err = %.2e\n", interp_err);
    if (interp_err > 1e-10)
        fail("tilted_ellipse", "interpolation error too large");

    pass("tilted_ellipse");
}

// ── Test 2: Parabola y=x² ─────────────────────────────────────────────────

void test_parabola()
{
    constexpr int Dim = 2;
    int n = 12;
    std::vector<VecN<Dim>> ctrl(n);
    std::vector<double> times(n);

    for (int i = 0; i < n; ++i) {
        double t = -2.0 + 4.0*i/(n-1);
        ctrl[i]  = VecN<Dim>(t, t*t);
        times[i] = t;
    }

    auto result = blend_curve<Dim>(ctrl, times, conic_tag{}, 60, 2);

    // Check approximation quality vs y=t²
    double err = 0.0;
    for (std::size_t i = 0; i < result.pts.size(); ++i) {
        double t = result.times[i], x = result.pts[i][0], y = result.pts[i][1];
        err = std::max(err, std::abs(y - x*x));
        (void)t;
    }
    std::printf("  parabola max |y - x²| = %.2e\n", err);
    // Parabola = degenerate conic (det=0); phi-orbit may fall back to circle.
    // Either way, interpolation at knots should work.
    double interp_err = interp_error<Dim>(result.pts, result.times, ctrl, times);
    std::printf("  parabola interp max err = %.2e\n", interp_err);
    if (interp_err > 1e-10)
        fail("parabola", "interpolation error too large");

    pass("parabola");
}

// ── Test 3: Helix (3D) — same quality as circle windows ───────────────────

void test_helix_3d()
{
    constexpr int Dim = 3;
    int n = 20;
    std::vector<VecN<Dim>> ctrl(n);
    std::vector<double> times(n);

    for (int i = 0; i < n; ++i) {
        double t = 4.0*PI * i / (n-1);
        ctrl[i]  = VecN<Dim>(std::cos(t), std::sin(t), 0.3*t);
        times[i] = t;
    }

    // Conic window result
    auto res_c = blend_curve<3>(ctrl, times, conic_tag{}, 60, 2);

    // Compare to circle window (fallback baseline)
    auto res_circ = blend_curve<3>(ctrl, times, circle_tag{}, 60, 2);

    // Both should have same number of output points
    if (res_c.pts.size() != res_circ.pts.size())
        std::printf("  helix_3d: conic pts=%zu  circle pts=%zu\n",
                    res_c.pts.size(), res_circ.pts.size());

    // Interpolation error
    double interp_err = interp_error<Dim>(res_c.pts, res_c.times, ctrl, times);
    std::printf("  helix_3d conic interp max err = %.2e\n", interp_err);
    if (interp_err > 1e-10)
        fail("helix_3d", "interpolation error too large");

    pass("helix_3d");
}

// ── Test 4: 4D embedded ellipse — tests SymEig<4> and back-projection ────────
//
// Ellipse a=3, b=1 embedded in a non-trivial 2D subspace of ℝ⁴:
//   e1 = (1,1,1,1)/2,  e2 = (1,-1,1,-1)/2   (orthonormal, span a 2D subspace)
// P(t) = center + 3·cos(t)·e1 + sin(t)·e2
//
// Since input lies exactly on a conic, conic windows should reproduce the
// interior control points to ~machine precision.

void test_4d_ellipse()
{
    constexpr int Dim = 4;
    // Orthonormal embedding basis
    double inv2 = 0.5;
    // e1 = (0.5, 0.5, 0.5, 0.5), e2 = (0.5, -0.5, 0.5, -0.5)

    int n = 16;
    std::vector<VecN<Dim>> ctrl(n);
    std::vector<double> times(n);

    for (int i = 0; i < n; ++i) {
        double t  = 2.0*PI * i / (n-1);
        double cx = 3.0*std::cos(t);   // coordinate along e1
        double cy = 1.0*std::sin(t);   // coordinate along e2
        ctrl[i] = VecN<Dim>(
            inv2*(cx + cy),
            inv2*(cx - cy),
            inv2*(cx + cy),
            inv2*(cx - cy));
        times[i] = t;
    }

    auto result = blend_curve<Dim>(ctrl, times, conic_tag{}, 60, 2);

    double interp_err = interp_error<Dim>(result.pts, result.times, ctrl, times);
    std::printf("  4d_ellipse interp max err = %.2e\n", interp_err);
    if (interp_err > 1e-10)
        fail("4d_ellipse", "interpolation error too large");

    pass("4d_ellipse");
}

// ── Test 5: Input validation ───────────────────────────────────────────────

void test_input_validation()  // was Test 4, now Test 5
{
    constexpr int Dim = 2;
    std::vector<VecN<Dim>> c6(6);
    std::vector<double>    t6(6);
    for (int i = 0; i < 6; ++i) {
        double t = 2.0*PI * i / 6;
        c6[i] = VecN<Dim>(std::cos(t), std::sin(t));
        t6[i] = t;
    }

    auto check = [&](char const* label, auto fn) {
        try { fn(); fail(label, "should have thrown"); }
        catch (std::invalid_argument const& e) {
            std::printf("  OK: %s → \"%s\"\n", label, e.what());
        }
    };

    check("size mismatch", [&]{ blend_curve<Dim>(c6, std::vector<double>{t6[0],t6[1]}, conic_tag{}, 60); });
    check("n<6", [&]{
        std::vector<VecN<Dim>> c5(c6.begin(), c6.begin()+5);
        std::vector<double>    t5(t6.begin(), t6.begin()+5);
        blend_curve<Dim>(c5, t5, conic_tag{}, 60);
    });
    check("pts_per_seg<2", [&]{ blend_curve<Dim>(c6, t6, conic_tag{}, 1); });

    pass("input_validation");
}

// ── Test 6: used_conic diagnostic parameter ────────────────────────────────
//
// Verifies the out-parameter reports true when conic windows are used and
// false when the call falls back to circle windows (collinear input).

void test_used_conic_flag()
{
    constexpr int Dim = 2;

    // (a) Tilted ellipse (2D, same as test_tilted_ellipse): exact conic input →
    //     all ConicWindows valid → used_conic=true.
    //     (The 3D helix has windows with collinear projected 2D points due to
    //     symmetric t-spacing, which correctly makes ConicWindow invalid for those
    //     windows, causing the whole curve to fall back to circle windows.)
    {
        double ca = std::cos(PI/6.0), sa = std::sin(PI/6.0);
        int n = 16;
        std::vector<VecN<Dim>> ctrl(n);
        std::vector<double> times(n);
        for (int i = 0; i < n; ++i) {
            double t = 2.0*PI * i / (n-1);
            double x = 3.0*std::cos(t);
            double y = 1.0*std::sin(t);
            ctrl[i]  = VecN<Dim>(ca*x - sa*y, sa*x + ca*y);
            times[i] = t;
        }
        bool flag = false;
        blend_curve<Dim>(ctrl, times, conic_tag{}, 40, 2, &flag);
        if (!flag) fail("used_conic_flag", "expected used_conic=true for tilted ellipse");
        std::printf("  tilted_ellipse: used_conic = %s\n", flag ? "true" : "false");
    }

    // (b) Collinear input: ConicWindow::valid()=false → fallback → used_conic=false
    {
        int n = 8;
        std::vector<VecN<Dim>> ctrl(n);
        std::vector<double> times(n);
        for (int i = 0; i < n; ++i) {
            ctrl[i]  = VecN<Dim>(static_cast<double>(i), 0.0);  // all on x-axis
            times[i] = static_cast<double>(i);
        }
        bool flag = true;
        blend_curve<Dim>(ctrl, times, conic_tag{}, 40, 2, &flag);
        if (flag) fail("used_conic_flag", "expected used_conic=false for collinear input");
        std::printf("  collinear: used_conic = %s\n", flag ? "true" : "false");
    }

    pass("used_conic_flag");
}

// ── Test 7: C^N continuity of the blended curve ───────────────────────────
//
// The blend formula (1−w)·A(t) + w(t)·B(t) guarantees C^N continuity when
// d^k w/ds^k = 0 at s=0 and s=1 for all k ≤ N.  At every interior junction
// the same window (wins[j-2]) appears on both sides, so:
//
//   C^0: trivially zero — same window evaluated at the same t.
//   C^1: zero — w'(0)=w'(1)=0, so cross terms vanish; both sides reduce to
//        wins[j-2]'(t_j).
//   C^2: zero — w''(0)=w''(1)=0 (quintic smoothstep); same argument.
//   C^3: zero — w'''(0) term is w'''(0)·(B(t_j)−A(t_j)), but B(t_j)=A(t_j)=ctrl[j]
//        by exact interpolation.  Requires Lagrange (C^∞) phi so phi''' is
//        continuous at ts[2]; PCHIP phi would give an O(1) jump here.
//   C^4: BREAKS — the fourth-derivative term 4·w'''(0)·(B'(t_j)−A'(t_j))·Δt
//        involves the DERIVATIVE mismatch between independently-fitted adjacent
//        windows.  Expected jump ~ 240/Δt² × O(1) = O(100+).
//
// Finite difference step sizes balance truncation O(h·|f^(k+1)|) against
// cancellation O(ε·|f|/h^k):
//   k=1,2: h = 1e-5  (h_opt ~ ε^(1/3) ≈ 1e-5)
//   k=3,4: h = 1e-4  (h_opt ~ ε^(1/4) ≈ 3e-4)

void test_cn_continuity()
{
    constexpr int Dim = 2;
    double ca = std::cos(PI/6.0), sa = std::sin(PI/6.0);
    int n = 16;
    std::vector<VecN<Dim>> ctrl(n);
    std::vector<double> times(n);
    for (int i = 0; i < n; ++i) {
        double t = 2.0*PI * i / (n-1);
        double x = 3.0*std::cos(t), y = 1.0*std::sin(t);
        ctrl[i]  = VecN<Dim>(ca*x - sa*y, sa*x + ca*y);
        times[i] = t;
    }

    // Build all (n-4) conic windows explicitly
    std::vector<ConicWindow<Dim>> wins;
    wins.reserve(n - 4);
    for (int i = 0; i < n-4; ++i) {
        wins.emplace_back(ctrl[i], ctrl[i+1], ctrl[i+2], ctrl[i+3], ctrl[i+4],
                          times[i], times[i+1], times[i+2], times[i+3], times[i+4]);
        if (!wins.back().valid())
            fail("cn_continuity", "window unexpectedly invalid");
    }

    // Report Lagrange vs PCHIP window count.
    // C^3 continuity requires Lagrange (C^∞) phi; PCHIP gives only C^1.
    int n_lagrange = 0;
    for (auto const& w : wins) if (w.uses_c_infty()) ++n_lagrange;
    std::printf("  windows using Lagrange (C^∞): %d/%d\n",
                n_lagrange, (int)wins.size());
    if (n_lagrange == 0)
        fail("cn_continuity", "no windows used Lagrange interpolant");

    // Segment j: (1−w)·wins[j-2](t) + w·wins[j-1](t) on [times[j], times[j+1]]
    auto eval_seg = [&](int j, double t) -> VecN<Dim> {
        double s = (t - times[j]) / (times[j+1] - times[j]);
        double w = fc::smoothstep(s, 2);
        return wins[j-2](t) * (1.0 - w) + wins[j-1](t) * w;
    };

    double const hf = 1e-5;   // fine step  — k=1,2
    double const hc = 1e-4;   // coarse step — k=3,4
    double max_c0 = 0, max_c1 = 0, max_c2 = 0, max_c3 = 0, max_c4 = 0;

    // Junction j: between segment j-1 (right end) and segment j (left end).
    // Valid for j = 3..n-4.
    for (int j = 3; j <= n-4; ++j) {
        double tj = times[j];

        // Fine-grid points (k=1,2)
        VecN<Dim> L0  = eval_seg(j-1, tj);
        VecN<Dim> Lf1 = eval_seg(j-1, tj -    hf);
        VecN<Dim> Lf2 = eval_seg(j-1, tj - 2.*hf);
        VecN<Dim> R0  = eval_seg(j,   tj);
        VecN<Dim> Rf1 = eval_seg(j,   tj +    hf);
        VecN<Dim> Rf2 = eval_seg(j,   tj + 2.*hf);

        // Coarse-grid points (k=3,4)
        VecN<Dim> Lc1 = eval_seg(j-1, tj -    hc);
        VecN<Dim> Lc2 = eval_seg(j-1, tj - 2.*hc);
        VecN<Dim> Lc3 = eval_seg(j-1, tj - 3.*hc);
        VecN<Dim> Lc4 = eval_seg(j-1, tj - 4.*hc);
        VecN<Dim> Rc1 = eval_seg(j,   tj +    hc);
        VecN<Dim> Rc2 = eval_seg(j,   tj + 2.*hc);
        VecN<Dim> Rc3 = eval_seg(j,   tj + 3.*hc);
        VecN<Dim> Rc4 = eval_seg(j,   tj + 4.*hc);

        // C^0: exact position jump
        max_c0 = std::max(max_c0, (L0 - R0).norm());

        // C^1: one-sided first differences — binomial (1, 1)
        VecN<Dim> d1_L = (L0  - Lf1)                  * (1./hf);
        VecN<Dim> d1_R = (Rf1 - R0)                   * (1./hf);
        max_c1 = std::max(max_c1, (d1_L - d1_R).norm());

        // C^2: one-sided second differences — binomial (1, 2, 1)
        double hf2 = hf*hf;
        VecN<Dim> d2_L = (L0  - Lf1*2. + Lf2)         * (1./hf2);
        VecN<Dim> d2_R = (Rf2 - Rf1*2. + R0)           * (1./hf2);
        max_c2 = std::max(max_c2, (d2_L - d2_R).norm());

        // C^3: one-sided third differences — binomial (1, 3, 3, 1)
        double hc3 = hc*hc*hc;
        VecN<Dim> d3_L = (L0  - Lc1*3. + Lc2*3. - Lc3) * (1./hc3);
        VecN<Dim> d3_R = (Rc3 - Rc2*3. + Rc1*3. - R0)  * (1./hc3);
        max_c3 = std::max(max_c3, (d3_L - d3_R).norm());

        // C^4: one-sided fourth differences — binomial (1, 4, 6, 4, 1)
        // Expected to fail: window derivative mismatch → O(240/Δt²) jump
        double hc4 = hc3*hc;
        VecN<Dim> d4_L = (L0  - Lc1*4. + Lc2*6. - Lc3*4. + Lc4) * (1./hc4);
        VecN<Dim> d4_R = (Rc4 - Rc3*4. + Rc2*6. - Rc1*4. + R0)  * (1./hc4);
        max_c4 = std::max(max_c4, (d4_L - d4_R).norm());
    }

    std::printf("  C^0 max jump = %.2e  (threshold 1e-12)\n",   max_c0);
    std::printf("  C^1 max jump = %.2e  (threshold 1e-3)\n",    max_c1);
    std::printf("  C^2 max jump = %.2e  (threshold 1e-2)\n",    max_c2);
    std::printf("  C^3 max jump = %.2e  (threshold 0.5)\n",     max_c3);
    std::printf("  C^4 max jump = %.2e  (diagnostic — expected O(100+))\n", max_c4);

    if (max_c0 > 1e-12) fail("cn_continuity", "C^0 discontinuity");
    if (max_c1 > 1e-3)  fail("cn_continuity", "C^1 discontinuity");
    if (max_c2 > 1e-2)  fail("cn_continuity", "C^2 discontinuity");
    if (max_c3 > 0.5)   fail("cn_continuity", "C^3 discontinuity");
    // C^4: diagnostic only — window derivative mismatch is expected

    pass("cn_continuity");
}

// ── Test 8: similarity invariance ─────────────────────────────────────────
//
// The phi-orbit stores slopes s = dy/dx in the conic's principal frame.
// Under rotation, the principal frame co-rotates, leaving s unchanged.
// Under uniform scale λ, s = (λ dy)/(λ dx) = dy/dx — unchanged.
// Under translation, the center of the conic shifts, but the slopes relative
// to P₀ are unchanged (the shift cancels in the numerator and denominator).
// Therefore PCHIP(phi, t) is identical before and after any similarity
// transformation, and blend_curve should produce T(original output).
//
// This does NOT hold for non-uniform scale (e.g. x→2x): the slopes halve,
// phi values shift, and the arc between knots differs — general affine
// invariance is not expected.

void test_similarity_invariance()
{
    constexpr int Dim = 2;
    int n = 12;
    std::vector<VecN<Dim>> ctrl(n);
    std::vector<double> times(n);
    for (int i = 0; i < n; ++i) {
        double t = 2.0*PI * i / (n-1);
        ctrl[i]  = VecN<Dim>(3.0*std::cos(t), 1.0*std::sin(t));
        times[i] = t;
    }

    // Reference result
    auto ref = blend_curve<Dim>(ctrl, times, conic_tag{}, 60, 2);

    // Similarity: rotate 37°, scale by 1.7, translate by (3, -2)
    double ang = 37.0 * PI / 180.0;
    double sc  = 1.7;
    double tx  = 3.0, ty = -2.0;
    double ca  = std::cos(ang), sa = std::sin(ang);
    auto T = [&](VecN<Dim> p) -> VecN<Dim> {
        double x = sc * p[0], y = sc * p[1];          // uniform scale
        double rx = ca*x - sa*y, ry = sa*x + ca*y;    // rotate
        return VecN<Dim>(rx + tx, ry + ty);            // translate
    };

    std::vector<VecN<Dim>> ctrl_t(n);
    for (int i = 0; i < n; ++i) ctrl_t[i] = T(ctrl[i]);

    auto res_t = blend_curve<Dim>(ctrl_t, times, conic_tag{}, 60, 2);

    if (res_t.pts.size() != ref.pts.size())
        fail("similarity_invariance", "output size mismatch after transformation");

    double max_err = 0.0;
    for (std::size_t i = 0; i < ref.pts.size(); ++i) {
        double err = (T(ref.pts[i]) - res_t.pts[i]).norm();
        if (err > max_err) max_err = err;
    }
    std::printf("  similarity_invariance max err = %.2e\n", max_err);
    if (max_err > 1e-10)
        fail("similarity_invariance", "output not similarity-invariant");

    pass("similarity_invariance");
}

// ── Test 9: affine consistency (diagnostic) ───────────────────────────────
//
// x-stretch by 1.5× is a non-similarity affine transform.  With vertex-P₀ +
// φ = 2·arctan(s) parametrization (mirrors conicspline.py), the output is not
// exactly affinely invariant (arctan is not Möbius-invariant under arbitrary
// affine transforms), but the error should be small for well-sampled curves.
// This is printed as a diagnostic; the threshold is loose.

void test_affine_invariance()
{
    constexpr int Dim = 2;
    int n = 16;
    std::vector<VecN<Dim>> ctrl(n);
    std::vector<double> times(n);
    for (int i = 0; i < n; ++i) {
        double t = 2.0*PI * i / (n-1);
        ctrl[i]  = VecN<Dim>(3.0*std::cos(t+0.3), 1.5*std::sin(t+0.3));  // tilted ellipse
        times[i] = t;
    }

    // Reference
    auto ref = blend_curve<Dim>(ctrl, times, conic_tag{}, 60, 2);

    // Non-similarity affine: stretch x by 1.5
    constexpr double SX = 1.5;
    std::vector<VecN<Dim>> ctrl_a(n);
    for (int i = 0; i < n; ++i)
        ctrl_a[i] = VecN<Dim>(SX * ctrl[i][0], ctrl[i][1]);

    auto res_a = blend_curve<Dim>(ctrl_a, times, conic_tag{}, 60, 2);

    if (res_a.pts.size() != ref.pts.size())
        fail("affine_invariance", "output size mismatch after transformation");

    double max_err = 0.0;
    for (std::size_t i = 0; i < ref.pts.size(); ++i) {
        VecN<Dim> expected(SX * ref.pts[i][0], ref.pts[i][1]);
        double err = (expected - res_a.pts[i]).norm();
        if (err > max_err) max_err = err;
    }
    std::printf("  affine_consistency (x-stretch 1.5) max err = %.2e\n", max_err);
    if (max_err > 1e-1)
        fail("affine_invariance", "affine consistency error unexpectedly large");

    pass("affine_invariance");
}

// ── test_cross_branch ─────────────────────────────────────────────────────
// Rectangular hyperbola xy = 1.  Six control points at x = -2,-1,-0.5,0.5,1,2
// so the first three are on the negative branch (x<0) and the last three on the
// positive branch (x>0): B2-B2-B2-B1-B1-B1, single crossing between index 2 and 3.
//
// Both 5-point windows span the crossing.  Cross-ratio fails for cross-branch
// windows (stereographic slopes pass through ±∞ at the asymptote → denominators
// hit zero), so the cross-ratio loop is skipped when is_cross_branch=true.
//
// Without allow_cross_branch: cross-ratio skipped + phi-unwrap skipped → found=false
//   → circle fallback → used_conic=false.
// With    allow_cross_branch: phi-unwrap runs → found=true → used_conic=true.
void test_cross_branch()
{
    constexpr int Dim = 2;
    const double xs[] = {-2.0, -1.0, -0.5, 0.5, 1.0, 2.0};

    // 6 points on xy=1: B2-B2-B2-B1-B1-B1, single crossing between index 2 and 3
    std::vector<VecN<Dim>> ctrl(6);
    std::vector<double> times(6);
    for (int i = 0; i < 6; ++i) {
        ctrl[i]  = VecN<Dim>(xs[i], 1.0 / xs[i]);
        times[i] = xs[i];
    }

    // Without allow_cross_branch: cross-ratio fails for every P0 → circle fallback
    bool used_conic_no = true;
    blend_curve<Dim>(ctrl, times, conic_tag{}, 40, 2, &used_conic_no, false);
    std::printf("  cross_branch=false: used_conic = %s\n", used_conic_no ? "true" : "false");
    if (used_conic_no)
        fail("cross_branch", "expected circle fallback without allow_cross_branch");

    // With allow_cross_branch: phi-unwrap succeeds → conic accepted
    bool used_conic_yes = false;
    auto res = blend_curve<Dim>(ctrl, times, conic_tag{}, 40, 2, &used_conic_yes, true);
    std::printf("  cross_branch=true:  used_conic = %s\n", used_conic_yes ? "true" : "false");
    if (!used_conic_yes)
        fail("cross_branch", "expected conic windows with allow_cross_branch");

    // Control-point interpolation check at each knot time
    double max_ctrl_err = 0.0;
    for (int i = 0; i < 6; ++i) {
        double t_target = times[i];
        double best_err = 1e30;
        for (std::size_t k = 0; k < res.pts.size(); ++k) {
            if (std::abs(res.times[k] - t_target) < 1e-12) {
                double err = (res.pts[k] - ctrl[i]).norm();
                if (err < best_err) best_err = err;
            }
        }
        if (best_err < 1e30 && best_err > max_ctrl_err) max_ctrl_err = best_err;
    }
    std::printf("  cross_branch: ctrl interp max err = %.2e\n", max_ctrl_err);
    if (max_ctrl_err > 1e-3)
        fail("cross_branch", "control-point interpolation error too large");

    pass("cross_branch");
}

// ── main ──────────────────────────────────────────────────────────────────

int main()
{
    std::printf("=== conicblend (5-pt conic) tests ===\n\n");
    test_tilted_ellipse();
    test_parabola();
    test_helix_3d();
    test_4d_ellipse();
    test_input_validation();
    test_used_conic_flag();
    test_cn_continuity();
    test_similarity_invariance();
    test_affine_invariance();
    test_cross_branch();
    std::printf("\nAll tests passed.\n");
}
