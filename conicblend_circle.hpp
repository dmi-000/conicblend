// conicblend_circle.hpp — C^N 3D/nD curve interpolation from circle arcs
// Header-only C++17 library.
//
// A circle is a special case of a conic (equal eigenvalues), so this header
// provides the simplest window type.  For the general 5-point conic-section
// windows see conicblend.hpp (not yet implemented).
//
// Architecture (mirrors conicspline's design):
//   • Every 3 consecutive control points define a unique circumscribed circle
//     in the plane they span — exact, no SVD projection needed.
//   • Each interior segment j is covered by two overlapping 3-point windows:
//       win[j-1]  through  ctrl[j-1], ctrl[j],   ctrl[j+1]
//       win[j]    through  ctrl[j],   ctrl[j+1],  ctrl[j+2]
//   • Blend:  (1-w)·win[j-1](t)  +  w·win[j](t),  w = smoothstep(s, N)
//   • C^N continuity: at t = times[j] the same window (win[j-1]) evaluates
//     from both sides → identical value and N vanishing derivatives.
//
// Torsion is implicit: adjacent windows live in slightly different planes,
// naturally encoding Frenet binormal rotation without any explicit τ estimate.
//
// Template parameter Dim selects the ambient dimension (≥ 2).
// The circumcenter formula is purely scalar (dot products); Gram-Schmidt
// replaces the cross product to complete the in-plane basis in any Dim.
//
// Minimum control points: 4  (yields 1 blended segment, j=1).
// Blended segments: j = 1, …, n−3.
//
// 3D convenience (no template syntax required):
//   using fc::Vec3;                          // = VecN<3>
//   auto r = fc::blend_curve(ctrl3d, times); // calls blend_curve<3>
//
// nD usage:
//   using V4 = fc::VecN<4>;
//   auto r = fc::blend_curve<4>(ctrl, times);
//
// Tag dispatch (forward-compatible with conicblend.hpp conic windows):
//   auto r = fc::blend_curve(ctrl, times, fc::circle_tag{});

#pragma once

// Define FC_NO_EXCEPTIONS before including this header to replace all
// std::invalid_argument throws with std::terminate().  Useful in
// -fno-exceptions builds or on embedded targets.
#ifndef FC_NO_EXCEPTIONS
#  include <stdexcept>
#  define FC_THROW(msg) throw std::invalid_argument(msg)
#else
#  include <cstdlib>
#  define FC_THROW(msg) std::terminate()
#endif

#include <array>
#include <cmath>
#include <type_traits>
#include <vector>

#if __cplusplus >= 202002L
#  include <numbers>
#  include <span>
#endif

namespace fc {

// ── Window-type tags for compile-time dispatch ────────────────────────────
//
//   fc::blend_curve(ctrl, times, fc::circle_tag{})  — this header
//   fc::blend_curve(ctrl, times, fc::conic_tag{})   — conicblend.hpp (future)
//
// Calling conic_tag{} without conicblend.hpp included gives a clear
// compile error rather than silently falling back to circle windows.
struct circle_tag {};
struct conic_tag  {};   // reserved; defined in conicblend.hpp

// ── Constants ──────────────────────────────────────────────────────────────

namespace detail {
#if __cplusplus >= 202002L
    constexpr double PI     = std::numbers::pi;
    constexpr double TWO_PI = 2.0 * std::numbers::pi;
#else
    constexpr double PI     = 3.14159265358979323846;
    constexpr double TWO_PI = 6.28318530717958647692;
#endif
}

// ── VecN<Dim>: fixed-size N-vector ────────────────────────────────────────

template <int Dim>
struct VecN {
    static_assert(Dim >= 2, "Dimension must be at least 2");

    std::array<double, Dim> data{};

    VecN() = default;

    // Variadic constructor: VecN<3>(x, y, z)  or  VecN<3>{x, y, z}
#if __cplusplus >= 202002L
    template<typename... Args> requires (sizeof...(Args) == static_cast<std::size_t>(Dim))
    VecN(Args... args) : data{static_cast<double>(args)...} {}
#else
    template<typename... Args,
             typename = std::enable_if_t<sizeof...(Args) == Dim>>
    VecN(Args... args) : data{static_cast<double>(args)...} {}
#endif

    double  operator[](int i) const { return data[i]; }
    double& operator[](int i)       { return data[i]; }

    VecN operator+(VecN const& b) const {
        VecN r; for (int i=0;i<Dim;++i) r[i]=data[i]+b[i]; return r;
    }
    VecN operator-(VecN const& b) const {
        VecN r; for (int i=0;i<Dim;++i) r[i]=data[i]-b[i]; return r;
    }
    VecN operator*(double s) const {
        VecN r; for (int i=0;i<Dim;++i) r[i]=data[i]*s; return r;
    }
    VecN& operator+=(VecN const& b) {
        for (int i=0;i<Dim;++i) data[i]+=b[i]; return *this;
    }

    double dot(VecN const& b) const {
        double s=0; for (int i=0;i<Dim;++i) s+=data[i]*b[i]; return s;
    }
    double norm2() const { return dot(*this); }
    double norm()  const { return std::sqrt(norm2()); }
    VecN normalized() const { return *this * (1.0/norm()); }
};

template <int Dim>
inline VecN<Dim> operator*(double s, VecN<Dim> const& v) { return v * s; }

// ── 3D convenience aliases ────────────────────────────────────────────────

using Vec3 = VecN<3>;

// ── CircleND: circumscribed circle of 3 points in ℝ^Dim ──────────────────

template <int Dim>
struct CircleND {
    VecN<Dim> center;
    double    radius;
    VecN<Dim> e1;   // unit in-plane vector along (p1−p0)
    VecN<Dim> e2;   // unit in-plane vector orthogonal to e1 (Gram-Schmidt)
};

using Circle3 = CircleND<3>;

// Fit the unique circumscribed circle through three non-collinear points.
//
// Circumcenter: C = p0 + s·a + t·b,  a = p1−p0,  b = p2−p0
//   [A  B][s]   [A/2]      A = a·a,  B = a·b,  C = b·b
//   [B  C][t] = [C/2]      D = AC − B²  (= 0 iff collinear, all dimensions)
//
// In-plane basis:
//   e1 = a / |a|
//   e2 = (b − (b·e1)·e1) / |...|   — Gram-Schmidt; replaces cross product in ℝ³
template <int Dim>
inline CircleND<Dim> circumcircle(
    VecN<Dim> const& p0,
    VecN<Dim> const& p1,
    VecN<Dim> const& p2)
{
    VecN<Dim> a = p1 - p0, b = p2 - p0;
    double A = a.dot(a), B = a.dot(b), C = b.dot(b);
    double D = A*C - B*B;

    if (D < 1e-14 * A * C)
        FC_THROW("fc::circumcircle: points are (near-)collinear; "
                 "add more control points or use a smaller window");

    double s = C*(A - B) / (2.0*D);
    double t = A*(C - B) / (2.0*D);

    VecN<Dim> center = p0 + s*a + t*b;
    double radius = (center - p0).norm();

    VecN<Dim> e1 = a.normalized();
    VecN<Dim> b_perp = b - b.dot(e1) * e1;
    VecN<Dim> e2 = b_perp.normalized();

#if __cplusplus >= 202002L
    return {.center = center, .radius = radius, .e1 = e1, .e2 = e2};
#else
    return {center, radius, e1, e2};
#endif
}

// ── CircleWindow: arc through 3 control points ────────────────────────────
//
// Normal mode: circumscribed circle arc, φ(t) quadratic Lagrange in angle.
//
// Near-collinear fallback (D < 1e-10·A·C, radius > ~10⁵ × chord):
//   Uses a direct quadratic Lagrange interpolant in position space.
//   This gives smooth, exactly-interpolating, nearly-linear behaviour at
//   straight sections of the curve.  Torsion is undefined there but also
//   negligible.  circumcircle() is never called, so no exception is thrown.
//
// The C^N continuity guarantee is unaffected: it depends only on the same
// window being evaluated from both sides of each knot, regardless of type.
template <int Dim>
class CircleWindow {
    bool          linear_;   // true = near-collinear, use Lagrange position
    VecN<Dim>     lp_[3];   // stored control points for linear mode
    CircleND<Dim> circ_;
    double t0_, t1_, t2_;
    double phi0_, phi1_, phi2_;

    static double wrap(double d) {
        return std::remainder(d, detail::TWO_PI);
    }

    double lagrange_basis(double t, int k) const {
        // L_k(t) for quadratic interpolation at t0_,t1_,t2_
        switch (k) {
            case 0: return (t-t1_)*(t-t2_) / ((t0_-t1_)*(t0_-t2_));
            case 1: return (t-t0_)*(t-t2_) / ((t1_-t0_)*(t1_-t2_));
            default:return (t-t0_)*(t-t1_) / ((t2_-t0_)*(t2_-t1_));
        }
    }

    double phi_of_t(double t) const {
        return phi0_*lagrange_basis(t,0)
             + phi1_*lagrange_basis(t,1)
             + phi2_*lagrange_basis(t,2);
    }

public:
    CircleWindow(
        VecN<Dim> const& p0, VecN<Dim> const& p1, VecN<Dim> const& p2,
        double t0, double t1, double t2)
        : linear_{false}, circ_{}, t0_{t0}, t1_{t1}, t2_{t2}
    {
        // Check near-collinearity before calling circumcircle.
        // D = AC − B² = |a×b|² (Lagrange identity); zero iff a ∥ b.
        // Threshold 1e-10·AC corresponds to sin(angle) < ~1e-5, i.e. the
        // circumradius > ~10⁵ × chord — genuinely straight for any practical curve.
        VecN<Dim> a = p1 - p0, b = p2 - p0;
        double A = a.dot(a), B = a.dot(b), C = b.dot(b);
        double D = A*C - B*B;
        if (D < 1e-10 * (A*C + 1e-30)) {
            linear_ = true;
            lp_[0] = p0;  lp_[1] = p1;  lp_[2] = p2;
            return;
        }

        circ_ = circumcircle<Dim>(p0, p1, p2);

        auto angle = [&](VecN<Dim> const& p) {
            VecN<Dim> r = p - circ_.center;
            return std::atan2(r.dot(circ_.e2), r.dot(circ_.e1));
        };

        phi0_ = angle(p0);
        double d1 = wrap(angle(p1) - phi0_);
        double d2 = wrap(angle(p2) - phi0_);

        // If d1 and d2 have opposite signs the arc crosses the ±π branch cut.
        if (d1 * d2 < 0.0)
            d2 += (d1 > 0.0 ? detail::TWO_PI : -detail::TWO_PI);

        phi1_ = phi0_ + d1;
        phi2_ = phi0_ + d2;
    }

    VecN<Dim> operator()(double t) const {
        if (linear_) {
            return lp_[0]*lagrange_basis(t,0)
                 + lp_[1]*lagrange_basis(t,1)
                 + lp_[2]*lagrange_basis(t,2);
        }
        double phi = phi_of_t(t);
        return circ_.center
             + (circ_.radius * std::cos(phi)) * circ_.e1
             + (circ_.radius * std::sin(phi)) * circ_.e2;
    }
};

// ── Smoothstep blend weight ───────────────────────────────────────────────
//
// w(s): degree-(2N+1) polynomial, w(0)=0, w(1)=1, w^(k)(0)=w^(k)(1)=0 for k=1..N.
// The blend (1−w)·A(t) + w·B(t) is C^N at s=0 and s=1.
inline double smoothstep(double s, int N = 2)
{
#if __cplusplus >= 202002L
    if (s <= 0.0) [[unlikely]] return 0.0;
    if (s >= 1.0) [[unlikely]] return 1.0;
#else
    if (s <= 0.0) return 0.0;
    if (s >= 1.0) return 1.0;
#endif
    switch (N) {
        case 1: return s*s*(3.0 - 2.0*s);
        case 2: return s*s*s*(10.0 + s*(-15.0 + 6.0*s));
        case 3: return s*s*s*s*(35.0 + s*(-84.0 + s*(70.0 - 20.0*s)));
        default:
            FC_THROW("fc::smoothstep: N must be 1, 2, or 3");
            return 0.0;
    }
}

// ── BlendResult / blend_curve ─────────────────────────────────────────────

template <int Dim>
struct BlendResultND {
    std::vector<VecN<Dim>> pts;
    std::vector<double>    times;
};

using BlendResult = BlendResultND<3>;   // 3D convenience alias

// Build and evaluate the blended curve.
// ctrl.size() == times.size() >= 4; times must be strictly increasing.
// Returns pts_per_seg points per interior segment; final point included once.
// The curve passes exactly through ctrl[1], …, ctrl[n-2].
template <int Dim>
inline BlendResultND<Dim> blend_curve(
#if __cplusplus >= 202002L
    std::span<VecN<Dim> const> ctrl,
    std::span<double const>    times,
#else
    std::vector<VecN<Dim>> const& ctrl,
    std::vector<double>    const& times,
#endif
    int pts_per_seg = 60,
    int smooth_N    = 2)
{
    int n = static_cast<int>(ctrl.size());
    if (n != static_cast<int>(times.size()))
        FC_THROW("fc::blend_curve: ctrl and times must have equal length");
    if (n < 4)
        FC_THROW("fc::blend_curve: need at least 4 control points");
    if (pts_per_seg < 2)
        FC_THROW("fc::blend_curve: pts_per_seg must be >= 2");

    // win[i] interpolates ctrl[i], ctrl[i+1], ctrl[i+2].
    std::vector<CircleWindow<Dim>> wins;
    wins.reserve(n - 2);
    for (int i = 0; i < n-2; ++i)
        wins.emplace_back(ctrl[i], ctrl[i+1], ctrl[i+2],
                          times[i], times[i+1], times[i+2]);

    // Segment j runs from times[j] to times[j+1], uses win[j-1] and win[j].
    // At s=0: w=0 → blend = win[j-1](times[j]) = ctrl[j]    ✓
    // At s=1: w=1 → blend = win[j]  (times[j+1]) = ctrl[j+1] ✓
    // C^N junction: both segments j-1 and j evaluate win[j-1] at times[j].
    int n_segs = n - 3;
    BlendResultND<Dim> out;
    out.pts.reserve(n_segs * pts_per_seg + 1);
    out.times.reserve(n_segs * pts_per_seg + 1);

    for (int j = 1; j <= n-3; ++j) {
        CircleWindow<Dim> const& wA = wins[j-1];
        CircleWindow<Dim> const& wB = wins[j];
        double t_lo = times[j], t_hi = times[j+1];

        int k_max = (j < n-3) ? pts_per_seg - 1 : pts_per_seg;
        for (int k = 0; k <= k_max; ++k) {
            double s = static_cast<double>(k) / pts_per_seg;
            double t = t_lo + s * (t_hi - t_lo);
            double w = smoothstep(s, smooth_N);
            VecN<Dim> A = wA(t), B = wB(t);
            out.pts.push_back(A*(1.0-w) + B*w);
            out.times.push_back(t);
        }
    }

    return out;
}

// ── Tagged overload: blend_curve(..., circle_tag{}) ───────────────────────

template <int Dim>
inline BlendResultND<Dim> blend_curve(
#if __cplusplus >= 202002L
    std::span<VecN<Dim> const> ctrl,
    std::span<double const>    times,
#else
    std::vector<VecN<Dim>> const& ctrl,
    std::vector<double>    const& times,
#endif
    circle_tag,
    int pts_per_seg = 60,
    int smooth_N    = 2)
{
    return blend_curve<Dim>(ctrl, times, pts_per_seg, smooth_N);
}

// ── 3D non-template overloads (no <3> required at call sites) ────────────

inline BlendResult blend_curve(
#if __cplusplus >= 202002L
    std::span<Vec3 const>   ctrl,
    std::span<double const> times,
#else
    std::vector<Vec3>   const& ctrl,
    std::vector<double> const& times,
#endif
    int pts_per_seg = 60,
    int smooth_N    = 2)
{
    return blend_curve<3>(ctrl, times, pts_per_seg, smooth_N);
}

inline BlendResult blend_curve(
#if __cplusplus >= 202002L
    std::span<Vec3 const>   ctrl,
    std::span<double const> times,
#else
    std::vector<Vec3>   const& ctrl,
    std::vector<double> const& times,
#endif
    circle_tag,
    int pts_per_seg = 60,
    int smooth_N    = 2)
{
    return blend_curve<3>(ctrl, times, pts_per_seg, smooth_N);
}

} // namespace fc
