// demo_scan.cpp — scan 2D test curves at varying n, output CSV for Python plotting.
//
// Compile: g++ -std=c++17 -O2 -I. -o demo_scan demo_scan.cpp
// Run:     ./demo_scan
//
// Outputs per curve {name}:
//   {name}_conv.csv      — n, err_conic, err_circle, pct_valid_wins
//   {name}_blend.csv     — t,x,y  (conic blend at N_DEMO ctrl pts)
//   {name}_circle.csv    — t,x,y  (circle blend at N_DEMO ctrl pts)
//   {name}_truth.csv     — t,x,y  (analytic, 500 pts)
//   {name}_ctrl.csv      — t,x,y  (N_DEMO control points)
//   {name}_wins.csv      — i,t_mid,valid  (per ConicWindow at N_DEMO)

#include "conicblend.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace fc;
using V2 = VecN<2>;

static constexpr double PI = 3.14159265358979323846;
static constexpr int    N_DEMO = 14;    // n used for visual demo plots

// ── curve definitions ──────────────────────────────────────────────────────

// Tilted ellipse: a=2, b=1, θ=30°
static V2 ellipse(double t) {
    double x =  2.0 * std::cos(t);
    double y =  1.0 * std::sin(t);
    double c = std::cos(PI / 6.0), s = std::sin(PI / 6.0);
    return V2(x*c - y*s, x*s + y*c);
}

// Parabola y = x²  (x = t)
static V2 parabola(double t) { return V2(t, t * t); }

// Lissajous 3:2
static V2 lissajous(double t) {
    return V2(2.0 * std::sin(2.0 * t), 2.0 * std::sin(3.0 * t + 0.3));
}

// Rose curve k=3
static V2 rose(double t) {
    double r = 2.0 * std::cos(3.0 * t);
    return V2(r * std::cos(t), r * std::sin(t));
}

// ── helpers ────────────────────────────────────────────────────────────────

static std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i)
        v[i] = a + (b - a) * i / (n - 1);
    return v;
}

using CurveFn = V2 (*)(double);

static std::vector<V2> sample_ctrl(CurveFn fn, double t0, double t1, int n) {
    auto ts = linspace(t0, t1, n);
    std::vector<V2> pts(n);
    for (int i = 0; i < n; ++i) pts[i] = fn(ts[i]);
    return pts;
}

static double max_err(BlendResultND<2> const& r, CurveFn truth) {
    double e = 0.0;
    for (std::size_t i = 0; i < r.pts.size(); ++i)
        e = std::max(e, (r.pts[i] - truth(r.times[i])).norm());
    return e;
}

// Count how many of the (n-4) ConicWindows are individually valid.
static int count_valid_wins(std::vector<V2> const& ctrl,
                            std::vector<double> const& ts) {
    int n = (int)ctrl.size(), cnt = 0;
    for (int i = 0; i < n - 4; ++i) {
        ConicWindow<2> w(ctrl[i], ctrl[i+1], ctrl[i+2], ctrl[i+3], ctrl[i+4],
                         ts[i], ts[i+1], ts[i+2], ts[i+3], ts[i+4]);
        if (w.valid()) ++cnt;
    }
    return cnt;
}

// ── CSV writers ────────────────────────────────────────────────────────────

static void write_pts2d(const char* fname, BlendResultND<2> const& r) {
    FILE* f = fopen(fname, "w");
    fprintf(f, "t,x,y\n");
    for (std::size_t i = 0; i < r.pts.size(); ++i)
        fprintf(f, "%.10f,%.10f,%.10f\n", r.times[i], r.pts[i][0], r.pts[i][1]);
    fclose(f);
}

static void write_truth(const char* fname, CurveFn fn, double t0, double t1) {
    FILE* f = fopen(fname, "w");
    fprintf(f, "t,x,y\n");
    for (int i = 0; i < 500; ++i) {
        double t = t0 + (t1 - t0) * i / 499.0;
        V2 p = fn(t);
        fprintf(f, "%.10f,%.10f,%.10f\n", t, p[0], p[1]);
    }
    fclose(f);
}

static void write_ctrl(const char* fname,
                       std::vector<V2> const& ctrl,
                       std::vector<double> const& ts) {
    FILE* f = fopen(fname, "w");
    fprintf(f, "t,x,y\n");
    for (std::size_t i = 0; i < ctrl.size(); ++i)
        fprintf(f, "%.10f,%.10f,%.10f\n", ts[i], ctrl[i][0], ctrl[i][1]);
    fclose(f);
}

static void write_windows(const char* fname,
                          std::vector<V2> const& ctrl,
                          std::vector<double> const& ts) {
    int n = (int)ctrl.size();
    FILE* f = fopen(fname, "w");
    fprintf(f, "i,t_mid,valid\n");
    for (int i = 0; i < n - 4; ++i) {
        ConicWindow<2> w(ctrl[i], ctrl[i+1], ctrl[i+2], ctrl[i+3], ctrl[i+4],
                         ts[i], ts[i+1], ts[i+2], ts[i+3], ts[i+4]);
        fprintf(f, "%d,%.10f,%d\n", i, ts[i+2], w.valid() ? 1 : 0);
    }
    fclose(f);
}

static void write_convergence(const char* fname, CurveFn fn, double t0, double t1) {
    static const int ns[] = {6, 8, 10, 12, 16, 20, 30, 40, 60};
    FILE* f = fopen(fname, "w");
    fprintf(f, "n,err_conic,err_circle,pct_valid_wins\n");
    for (int n : ns) {
        auto ts   = linspace(t0, t1, n);
        auto ctrl = sample_ctrl(fn, t0, t1, n);

        bool used_conic = false;
        auto rc = blend_curve<2>(ctrl, ts, conic_tag{},   40, 2, &used_conic);
        auto rk = blend_curve<2>(ctrl, ts, circle_tag{},  40, 2);

        double ec = max_err(rc, fn);
        double ek = max_err(rk, fn);
        int valid = count_valid_wins(ctrl, ts);
        double pct = 100.0 * valid / (n - 4);

        fprintf(f, "%d,%.6e,%.6e,%.1f\n", n, ec, ek, pct);
    }
    fclose(f);
}

// ── main ───────────────────────────────────────────────────────────────────

struct Curve { const char* name; CurveFn fn; double t0, t1; };

int main() {
    const Curve curves[] = {
        { "ellipse",   ellipse,   0.0,        2.0*PI     },
        { "parabola",  parabola,  -2.0,        2.0        },
        { "lissajous", lissajous,  0.0,        1.8*PI     },
        { "rose",      rose,       0.1,        PI - 0.1   },
    };

    printf("%-12s  %3s  %10s  %10s  %8s\n",
           "curve", "n", "err_conic", "err_circle", "conic_wins");
    printf("%s\n", std::string(52, '-').c_str());

    char buf[128];
    for (auto const& c : curves) {
        // Convergence
        snprintf(buf, sizeof(buf), "%s_conv.csv",   c.name);  write_convergence(buf, c.fn, c.t0, c.t1);
        snprintf(buf, sizeof(buf), "%s_truth.csv",  c.name);  write_truth(buf, c.fn, c.t0, c.t1);

        // Demo at N_DEMO
        auto ts   = linspace(c.t0, c.t1, N_DEMO);
        auto ctrl = sample_ctrl(c.fn, c.t0, c.t1, N_DEMO);

        bool used_conic = false;
        auto rc = blend_curve<2>(ctrl, ts, conic_tag{},  40, 2, &used_conic);
        auto rk = blend_curve<2>(ctrl, ts, circle_tag{}, 40, 2);

        snprintf(buf, sizeof(buf), "%s_blend.csv",  c.name);  write_pts2d(buf, rc);
        snprintf(buf, sizeof(buf), "%s_circle.csv", c.name);  write_pts2d(buf, rk);
        snprintf(buf, sizeof(buf), "%s_ctrl.csv",   c.name);  write_ctrl(buf, ctrl, ts);
        snprintf(buf, sizeof(buf), "%s_wins.csv",   c.name);  write_windows(buf, ctrl, ts);

        double ec = max_err(rc, c.fn);
        double ek = max_err(rk, c.fn);
        int    valid = count_valid_wins(ctrl, ts);
        printf("%-12s  %3d  %10.3e  %10.3e  %d/%d\n",
               c.name, N_DEMO, ec, ek, valid, N_DEMO - 4);
    }
    printf("\nAll CSV files written.\n");
}
