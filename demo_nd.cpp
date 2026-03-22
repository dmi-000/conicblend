// demo_nd.cpp — tests for conicblend_circle.hpp (nD and tag-dispatch paths)
//
// Compile:  clang++ -std=c++17 -O2 -I. -o demo_nd demo_nd.cpp
//
// Tests
//   1. Dim=2, unit circle:  exact reproduction to machine precision (~1e-15)
//   2. Dim=3, helix:        max error vs analytic helix; exact interpolation ~1e-15
//   3. Dim=4, torus curve:  exact interpolation at control points (~1e-15)
//   4. Collinearity guard:  circumcircle throws for Dim=2 and Dim=4
//   5. Input validation:    size mismatch, n<4, pts_per_seg<2, bad N

#include "conicblend_circle.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>

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

static void pass(char const* name) { std::printf("PASS  %s\n", name); }
static void fail(char const* name, char const* msg) {
    std::printf("FAIL  %s: %s\n", name, msg);
    std::exit(1);
}

// ── Test 1: Dim=2, unit circle ─────────────────────────────────────────────

void test_circle_2d()
{
    constexpr int Dim = 2;
    int n = 10;
    std::vector<VecN<Dim>> ctrl(n);
    std::vector<double>    times(n);
    for (int i = 0; i < n; ++i) {
        double t = 2.0*PI * i / n;
        ctrl[i]  = VecN<Dim>(std::cos(t), std::sin(t));
        times[i] = t;
    }

    auto result = blend_curve<Dim>(ctrl, times, 80, 2);

    auto truth = [](double t) -> VecN<Dim> {
        return VecN<Dim>(std::cos(t), std::sin(t));
    };
    double err = max_error<Dim>(result.pts, result.times, truth);
    std::printf("  Dim=2 circle max err = %.2e\n", err);
    if (err > 1e-12) fail("circle_2d", "error too large");

    // Tagged overload must produce identical output.
    auto r_tagged = blend_curve<Dim>(ctrl, times, circle_tag{}, 80, 2);
    assert(r_tagged.pts.size() == result.pts.size());
    for (std::size_t i = 0; i < result.pts.size(); ++i)
        assert((result.pts[i] - r_tagged.pts[i]).norm() == 0.0);

    pass("circle_2d");
}

// ── Test 2: Dim=3, helix — error vs analytic + exact interpolation ─────────
//
// Helix: P(t) = (cos t, sin t, 0.3 t).
// Tests two things:
//   a. Max position error vs analytic helix (approximation quality).
//   b. Exact interpolation: blended curve passes through every control point
//      to floating-point precision (algebraic identity, not approximation).

void test_helix_3d()
{
    constexpr int Dim = 3;
    int n = 20;
    std::vector<VecN<Dim>> ctrl(n);
    std::vector<double>    times(n);

    for (int i = 0; i < n; ++i) {
        double t = 4.0*PI * i / (n-1);
        ctrl[i]  = VecN<Dim>(std::cos(t), std::sin(t), 0.3*t);
        times[i] = t;
    }

    auto result = blend_curve<3>(ctrl, times, 60, 2);

    // (a) approximation error vs analytic
    auto truth = [](double t) -> VecN<3> {
        return VecN<3>(std::cos(t), std::sin(t), 0.3*t);
    };
    double approx_err = max_error<3>(result.pts, result.times, truth);
    std::printf("  Dim=3 helix approx max err = %.2e\n", approx_err);
    if (approx_err > 0.1) fail("helix_3d", "approximation error too large");

    // (b) exact interpolation: find each knot in the output and check distance
    double max_interp_err = 0.0;
    for (int j = 1; j <= n-2; ++j) {
        double tj = times[j];
        double best_d = 1e30; std::size_t best_i = 0;
        for (std::size_t i = 0; i < result.times.size(); ++i) {
            double d = std::abs(result.times[i] - tj);
            if (d < best_d) { best_d = d; best_i = i; }
        }
        double dt_step = (times[n-1]-times[0]) / ((n-3)*60.0);
        if (best_d > 2.0*dt_step) continue;
        max_interp_err = std::max(max_interp_err,
                                  (result.pts[best_i] - ctrl[j]).norm());
    }
    std::printf("  Dim=3 helix interp max err = %.2e\n", max_interp_err);
    if (max_interp_err > 1e-12) fail("helix_3d", "interpolation error too large");

    // 3D non-template overload must match the template call exactly.
    std::vector<Vec3> ctrl3(n);
    for (int i = 0; i < n; ++i)
        ctrl3[i] = Vec3(ctrl[i][0], ctrl[i][1], ctrl[i][2]);
    auto r3 = blend_curve(ctrl3, times, 60, 2);   // non-template 3D overload
    assert(r3.pts.size() == result.pts.size());
    for (std::size_t i = 0; i < result.pts.size(); ++i) {
        double d = (r3.pts[i] - result.pts[i]).norm();
        assert(d == 0.0);
    }

    pass("helix_3d");
}

// ── Test 3: Dim=4, 4D torus curve — exact interpolation ───────────────────

void test_4d_curve()
{
    constexpr int Dim = 4;
    int n = 16;
    std::vector<VecN<Dim>> ctrl(n);
    std::vector<double>    times(n);

    for (int i = 0; i < n; ++i) {
        double t  = 2.0*PI * i / (n-1);
        ctrl[i]   = VecN<Dim>(std::cos(t), std::sin(t),
                                std::cos(2*t), std::sin(2*t));
        times[i]  = t;
    }

    auto result = blend_curve<Dim>(ctrl, times, 60, 2);

    double max_interp_err = 0.0;
    for (int j = 1; j <= n-2; ++j) {
        double tj = times[j];
        double best_d = 1e30; std::size_t best_i = 0;
        for (std::size_t i = 0; i < result.times.size(); ++i) {
            double d = std::abs(result.times[i] - tj);
            if (d < best_d) { best_d = d; best_i = i; }
        }
        double dt_step = (times[n-1]-times[0]) / ((n-3)*60.0);
        if (best_d > 2.0*dt_step) continue;
        max_interp_err = std::max(max_interp_err,
                                  (result.pts[best_i] - ctrl[j]).norm());
    }
    std::printf("  Dim=4 exact interp max err = %.2e\n", max_interp_err);
    if (max_interp_err > 1e-12) fail("4d_curve", "interpolation error too large");
    pass("4d_curve");
}

// ── Test 4: collinearity guard ─────────────────────────────────────────────

void test_collinearity_guard()
{
    try {
        VecN<2> p0(0,0), p1(1,0), p2(2,0);
        circumcircle<2>(p0, p1, p2);
        fail("collinearity_2d", "should have thrown");
    } catch (std::invalid_argument const& e) {
        std::printf("  Dim=2 collinear exception: %s\n", e.what());
    }

    try {
        VecN<4> p0(0,0,0,0), p1(1,0,0,0), p2(2,0,0,0);
        circumcircle<4>(p0, p1, p2);
        fail("collinearity_4d", "should have thrown");
    } catch (std::invalid_argument const& e) {
        std::printf("  Dim=4 collinear exception: %s\n", e.what());
    }

    pass("collinearity_guard");
}

// ── Test 5: near-collinear windows — linear fallback, no exception ─────────
//
// A curve with a straight section: circle → straight line → circle.
// Near-collinear 3-point windows should silently use the Lagrange interpolant.
// The curve must (a) not throw, (b) pass through interior control points,
// (c) produce sensible (non-NaN) output.

void test_near_collinear()
{
    constexpr int Dim = 2;
    // 8 control points: curve bends, goes nearly straight, bends again.
    // Middle 3 points are nearly collinear (on y=0, x=1,2,3).
    std::vector<VecN<Dim>> ctrl = {
        VecN<Dim>(0.0,  1.0),   // 0: bend
        VecN<Dim>(0.5,  0.1),   // 1
        VecN<Dim>(1.0,  0.0),   // 2: nearly straight starts
        VecN<Dim>(2.0,  0.0),   // 3: exactly on x-axis
        VecN<Dim>(3.0,  0.0),   // 4: exactly on x-axis
        VecN<Dim>(3.5,  0.1),   // 5
        VecN<Dim>(4.0,  1.0),   // 6: bend
        VecN<Dim>(4.5,  2.0),   // 7
    };
    std::vector<double> times = {0,1,2,3,4,5,6,7};

    // Must not throw
    BlendResultND<Dim> result;
    try {
        result = blend_curve<Dim>(ctrl, times, 40, 2);
    } catch (std::exception const& e) {
        fail("near_collinear", e.what());
    }

    // No NaN in output
    for (std::size_t i = 0; i < result.pts.size(); ++i) {
        if (!std::isfinite(result.pts[i][0]) || !std::isfinite(result.pts[i][1]))
            fail("near_collinear", "NaN or Inf in output");
    }

    // Passes through interior control points (j=1..n-2 = j=1..6)
    int n = static_cast<int>(ctrl.size());
    double max_interp_err = 0.0;
    for (int j = 1; j <= n-2; ++j) {
        double tj = times[j];
        std::size_t best_i = 0;
        double best_d = 1e30;
        for (std::size_t i = 0; i < result.times.size(); ++i) {
            double d = std::abs(result.times[i] - tj);
            if (d < best_d) { best_d = d; best_i = i; }
        }
        double err = (result.pts[best_i] - ctrl[j]).norm();
        if (err > max_interp_err) max_interp_err = err;
    }
    std::printf("  near_collinear interp max err = %.2e\n", max_interp_err);
    if (max_interp_err > 1e-10)
        fail("near_collinear", "interpolation error too large");

    pass("near_collinear");
}

// ── Test 6: input validation ───────────────────────────────────────────────

void test_input_validation()  // was Test 5, now Test 6
{
    constexpr int Dim = 2;
    std::vector<VecN<Dim>> c5(5);
    std::vector<double>    t5(5);
    for (int i = 0; i < 5; ++i) {
        double t = 2.0*PI * i / 5;
        c5[i] = VecN<Dim>(std::cos(t), std::sin(t));
        t5[i] = t;
    }

    auto check = [&](char const* label, auto fn) {
        try { fn(); fail(label, "should have thrown"); }
        catch (std::invalid_argument const& e) {
            std::printf("  OK: %s → \"%s\"\n", label, e.what());
        }
    };

    check("size mismatch", [&]{ blend_curve<Dim>(c5, std::vector<double>{t5[0],t5[1],t5[2]}, 60); });
    check("n<4",           [&]{
        std::vector<VecN<Dim>> c3(c5.begin(), c5.begin()+3);
        std::vector<double>    t3(t5.begin(), t5.begin()+3);
        blend_curve<Dim>(c3, t3, 60);
    });
    check("pts_per_seg<2", [&]{ blend_curve<Dim>(c5, t5, 1); });
    check("bad N",         [&]{ blend_curve<Dim>(c5, t5, 60, 4); });

    pass("input_validation");
}

// ── main ──────────────────────────────────────────────────────────────────

int main()
{
    std::printf("=== conicblend_circle nD tests ===\n\n");
    test_circle_2d();
    test_helix_3d();
    test_4d_curve();
    test_collinearity_guard();
    test_near_collinear();
    test_input_validation();
    std::printf("\nAll tests passed.\n");
}
