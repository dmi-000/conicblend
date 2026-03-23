// demo_cylinder.cpp — compare ConicWindow vs CylinderWindow on 3D curves
//
// Test A: Helix r(t) = (cos 2t, sin 2t, 0.5t)
//   Lies exactly on a cylinder.  CylinderWindow should give near-zero knot
//   error; ConicWindow introduces O(κτh³) torsion error.
//
// Test B: Twisted cubic r(t) = (t, t², t³)
//   Not on any cylinder but has non-zero torsion.
//
// For each curve and each method: reports max deviation from the true curve
// at 200 intermediate sample points per segment.
//
// Compile:
//   g++ -std=c++17 -O2 -I. -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 \
//       -o demo_cylinder demo_cylinder.cpp

#include "conicblend_cylinder.hpp"
#include <cstdio>
#include <cmath>
#include <vector>

using namespace fc;
using V3 = VecN<3>;

// ── True curves ────────────────────────────────────────────────────────────
static V3 helix   (double t) { return V3{std::cos(2*t), std::sin(2*t), 0.5*t}; }
static V3 tcubic  (double t) { return V3{t, t*t, t*t*t}; }

// ── Max deviation from true curve at N interior samples ────────────────────
template <typename TrueFn>
static double max_dev(const BlendResult& r, TrueFn fn)
{
    double mx = 0.0;
    for (size_t i = 0; i < r.pts.size(); ++i) {
        V3 err = r.pts[i] - fn(r.times[i]);
        mx = std::max(mx, err.norm());
    }
    return mx;
}

// ── Run one curve ──────────────────────────────────────────────────────────
template <typename TrueFn>
static void run(const char* label, TrueFn fn,
                double t0, double t1, int n)
{
    std::printf("%-22s  n=%2d\n", label, n);

    std::vector<V3>     ctrl(n);
    std::vector<double> times(n);
    for (int i = 0; i < n; ++i) {
        times[i] = t0 + (t1-t0)*(double)i/(n-1);
        ctrl[i]  = fn(times[i]);
    }

    // ConicWindow (plane projection)
    bool used_c = false;
    auto r_conic = blend_curve(ctrl, times, conic_tag{}, 200, 2, &used_c);
    double dev_conic = max_dev(r_conic, fn);

    // CylinderWindow
    bool used_cyl = false;
    auto r_cyl = blend_curve(ctrl, times, cylinder_tag{}, 200, 2, &used_cyl);
    double dev_cyl = max_dev(r_cyl, fn);

    std::printf("  conic    (used=%d)  max_dev = %.3e\n", used_c,   dev_conic);
    std::printf("  cylinder (used=%d)  max_dev = %.3e\n", used_cyl, dev_cyl);
    if (dev_conic > 1e-15)
        std::printf("  improvement: %.1fx\n", dev_conic / (dev_cyl + 1e-20));
    std::printf("\n");
}

int main()
{
    std::printf("=== Test A: Helix r(t)=(cos 2t, sin 2t, 0.5t) ===\n\n");
    run("helix n=8",  helix, 0.0, 2*M_PI, 8);
    run("helix n=12", helix, 0.0, 2*M_PI, 12);
    run("helix n=20", helix, 0.0, 2*M_PI, 20);

    std::printf("=== Test B: Twisted cubic r(t)=(t,t²,t³) ===\n\n");
    run("tcubic n=8",  tcubic, 0.0, 1.0, 8);
    run("tcubic n=12", tcubic, 0.0, 1.0, 12);
    run("tcubic n=20", tcubic, 0.0, 1.0, 20);

    return 0;
}
