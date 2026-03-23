// conicblend.hpp — C^N curve interpolation using 5-point conic-section windows
// Header-only C++17 library.
//
// Architecture (mirrors conicspline's design):
//   • Every 5 consecutive control points define a best-fit conic arc in their
//     plane (SVD null-space of the 5×6 design matrix).
//   • Each interior segment j is covered by two overlapping 5-point windows:
//       win[j-2]  through  ctrl[j-2..j+2]
//       win[j-1]  through  ctrl[j-1..j+3]
//   • Blend:  (1-w)·win[j-2](t)  +  w·win[j-1](t),  w = smoothstep(s, N)
//   • C^N continuity: same window (win[j-2]) evaluates from both sides of t=times[j].
//
// Parametrization: φ = 2·arctan(s) orbit (mirrors conicspline.py) —
//   P₀ is the algebraic vertex of the conic (where ∂Q/∂x = 0), or the inner-3
//   control point minimising |L/M| as fallback.  Slopes s_i = (y_i−y₀)/(x_i−x₀)
//   give φ_i = 2·arctan(s_i); after ±π unwrapping, φ is monotone for valid windows.
//   A degree-4 Lagrange polynomial (or PCHIP fallback) interpolates φ(t); the
//   rational conic map x=x₀+u(s), y=y₀+s·u(s) evaluates the orbit.
//
// Using the algebraic vertex as P₀ ensures the same P₀ for any window on the
// same conic, so overlapping windows share the same orbit → exact blend for
// locally-conic data.  Works uniformly for ellipses, hyperbolas, and parabolas.
//
// Fallback: ConicWindows that fail (non-monotone φ, orbit reconstruction
//   error > 1e-3, or cross-branch without caller permission) fall back to
//   LagrangeWindow (degree-4 Lagrange through 5 control points).
//
// Minimum control points: 6  (yields 1 blended segment, j=2).
// Blended segments: j = 2, …, n−4.
//
// This file requires conicblend_circle.hpp (included automatically below).
//
// Usage:
//   auto r = fc::blend_curve<3>(ctrl, times, fc::conic_tag{});

#pragma once
#include "conicblend_circle.hpp"

#include <algorithm>
#include <cstring>
#include <variant>


namespace fc {

// ── Template Jacobi eigendecomposition of an N×N symmetric matrix ─────────
//
// Cyclic Jacobi sweeps: repeatedly zero the off-diagonal entry (p,q) with
// largest absolute value.  Converges to machine precision in ~N² sweeps.
//
// Returns eigenvalues sorted ascending in vals[], eigenvectors as columns
// of vecs[N][N] (vecs[col][row] indexing for column-major access).
//
// Used for:
//   N=6: smallest eigenvector of M^T M  →  conic coefficients [A..F]
//   N=Dim: two largest eigenvectors of covariance  →  best-fit 2D plane

namespace detail {

template <int N>
struct SymEig {
    double val[N];
    double vec[N][N];  // vec[j] = j-th eigenvector (column j)

    // Compute eigendecomposition of symmetric matrix A (row-major, A[row][col]).
    // A is modified in-place (diagonalized).
    void compute(double A[N][N])
    {
        // Initialize eigenvectors to identity
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                vec[i][j] = (i == j) ? 1.0 : 0.0;

        for (int sweep = 0; sweep < 50 * N; ++sweep) {
            // Find largest off-diagonal element
            double off = 0.0;
            for (int i = 0; i < N-1; ++i)
                for (int j = i+1; j < N; ++j)
                    off += A[i][j] * A[i][j];
            if (off < 1e-28) break;

            // Cyclic sweeps over all (p,q) pairs
            for (int p = 0; p < N-1; ++p) {
                for (int q = p+1; q < N; ++q) {
                    if (std::abs(A[p][q]) < 1e-30) continue;
                    double diff = A[q][q] - A[p][p];
                    double t_rot, c, s;
                    if (std::abs(A[p][q]) < 1e-15 * std::abs(diff)) {
                        t_rot = A[p][q] / diff;
                    } else {
                        double theta = 0.5 * diff / A[p][q];
                        t_rot = 1.0 / (std::abs(theta) + std::sqrt(1.0 + theta*theta));
                        if (theta < 0.0) t_rot = -t_rot;
                    }
                    c = 1.0 / std::sqrt(1.0 + t_rot*t_rot);
                    s = t_rot * c;
                    double tau = s / (1.0 + c);
                    double h = t_rot * A[p][q];
                    A[p][p] -= h;
                    A[q][q] += h;
                    A[p][q] = 0.0;  A[q][p] = 0.0;
                    for (int r = 0; r < N; ++r) {
                        if (r != p && r != q) {
                            double apq_r = A[p][r], aqr = A[q][r];
                            A[p][r] = apq_r - s*(aqr + tau*apq_r);
                            A[q][r] = aqr   + s*(apq_r - tau*aqr);
                            A[r][p] = A[p][r];
                            A[r][q] = A[q][r];
                        }
                        // Update eigenvector matrix (columns rotate)
                        double vp = vec[p][r], vq = vec[q][r];
                        vec[p][r] = vp - s*(vq + tau*vp);
                        vec[q][r] = vq + s*(vp - tau*vq);
                    }
                }
            }
        }
        // Copy diagonal to val
        for (int i = 0; i < N; ++i) val[i] = A[i][i];

        // Sort eigenvalues ascending, permute eigenvectors accordingly
        for (int i = 0; i < N-1; ++i) {
            int mn = i;
            for (int j = i+1; j < N; ++j)
                if (val[j] < val[mn]) mn = j;
            if (mn != i) {
                std::swap(val[i], val[mn]);
                for (int r = 0; r < N; ++r) std::swap(vec[i][r], vec[mn][r]);
            }
        }
    }
};

// ── Fit conic to exactly 5 2D points ─────────────────────────────────────
//
// Design matrix M[i] = [x²,x·y,y²,x,y,1] (5×6).
// We want the null vector of M — the smallest eigenvector of M^T·M (6×6).
// Returns [A,B,C,D,E,F] in coeffs[6].  Not normalized.

inline void fit_conic_2d(double const pts[5][2], double coeffs[6])
{
    // Build M^T * M  (6×6 symmetric)
    double MtM[6][6] = {};
    for (int i = 0; i < 5; ++i) {
        double x = pts[i][0], y = pts[i][1];
        double row[6] = {x*x, x*y, y*y, x, y, 1.0};
        for (int r = 0; r < 6; ++r)
            for (int c = 0; c < 6; ++c)
                MtM[r][c] += row[r] * row[c];
    }
    SymEig<6> eig;
    eig.compute(MtM);
    // Smallest eigenvalue → conic null vector
    for (int i = 0; i < 6; ++i) coeffs[i] = eig.vec[0][i];
}

// ── Closed-form 2×2 symmetric eigensolver ────────────────────────────────
//
// Matrix [[a, b], [b, c]].
// Eigenvalues: (a+c)/2 ± sqrt(((a-c)/2)^2 + b^2).
// Eigenvectors computed from (A - λ·I)v = 0.

struct Eig2x2 {
    double val[2];     // ascending
    double vec[2][2];  // vec[j] = j-th eigenvector
};

inline Eig2x2 eigh2(double a, double b, double c)
{
    Eig2x2 r;
    double mid = 0.5*(a + c);
    double d   = std::sqrt(0.25*(a-c)*(a-c) + b*b);
    r.val[0] = mid - d;  // smaller
    r.val[1] = mid + d;

    // Eigenvectors: for each eigenvalue λ, solve (a-λ)v₀ + b·v₁ = 0
    for (int k = 0; k < 2; ++k) {
        double lam = r.val[k];
        double v0 = b, v1 = lam - a;  // from second row: b·v₀ + (c-λ)v₁=0 → scaled
        // If b≈0 use first row: (a-λ)v₀=0 → v₀=0 or eigenvalue
        double norm = std::sqrt(v0*v0 + v1*v1);
        if (norm < 1e-15) {
            // Degenerate (a=c=λ, b=0): identity direction
            r.vec[k][0] = (k == 0) ? 1.0 : 0.0;
            r.vec[k][1] = (k == 0) ? 0.0 : 1.0;
        } else {
            r.vec[k][0] = v0 / norm;
            r.vec[k][1] = v1 / norm;
        }
    }
    return r;
}

// ── Fixed-size PCHIP interpolant for n=5 nodes ───────────────────────────
//
// Fritsch-Carlson monotone cubic Hermite interpolation.
// Stores nodes (t[5], y[5]) and computed tangents d[5].
// Evaluates pos, first derivative, and second derivative.

struct Pchip5 {
    double t[5], y[5], d[5];

    // Initialize from 5 (time, value) pairs.
    // t[] must be strictly increasing.
    void init(double const ts[5], double const ys[5])
    {
        for (int i = 0; i < 5; ++i) { t[i] = ts[i]; y[i] = ys[i]; }

        double h[4], delta[4];
        for (int k = 0; k < 4; ++k) {
            h[k]     = t[k+1] - t[k];
            delta[k] = (y[k+1] - y[k]) / h[k];
        }

        // Initialize tangents (three-point finite difference at interior,
        // one-sided at endpoints)
        d[0] = delta[0];
        d[4] = delta[3];
        for (int k = 1; k <= 3; ++k) {
            if (delta[k-1] * delta[k] <= 0.0)
                d[k] = 0.0;
            else
                d[k] = (delta[k-1] + delta[k]) * 0.5;
        }

        // Fritsch-Carlson monotonicity correction
        for (int k = 0; k < 4; ++k) {
            if (std::abs(delta[k]) < 1e-30) {
                d[k] = 0.0;  d[k+1] = 0.0;
                continue;
            }
            double alpha = d[k]   / delta[k];
            double beta  = d[k+1] / delta[k];
            double rho   = alpha*alpha + beta*beta;
            if (rho > 9.0) {
                double tau = 3.0 / std::sqrt(rho);
                d[k]   = tau * alpha * delta[k];
                d[k+1] = tau * beta  * delta[k];
            }
        }
    }

    // Find interval k such that t[k] ≤ s ≤ t[k+1], returns local τ ∈ [0,1].
    int locate(double s, double& tau) const
    {
        int k = 0;
        for (int i = 0; i < 3; ++i)
            if (s >= t[i+1]) k = i+1;
        double hk = t[k+1] - t[k];
        tau = (hk > 1e-30) ? (s - t[k]) / hk : 0.0;
        tau = std::max(0.0, std::min(1.0, tau));
        return k;
    }

    double eval(double s) const
    {
        double tau; int k = locate(s, tau);
        double tau2 = tau*tau, tau3 = tau2*tau;
        double h00 = 2*tau3 - 3*tau2 + 1;
        double h10 = tau3 - 2*tau2 + tau;
        double h01 = -2*tau3 + 3*tau2;
        double h11 = tau3 - tau2;
        double hk = t[k+1] - t[k];
        return y[k]*h00 + hk*d[k]*h10 + y[k+1]*h01 + hk*d[k+1]*h11;
    }

    // First derivative dy/dt
    double deriv1(double s) const
    {
        double tau; int k = locate(s, tau);
        double tau2 = tau*tau;
        double dh00 = 6*tau2 - 6*tau;
        double dh10 = 3*tau2 - 4*tau + 1;
        double dh01 = -6*tau2 + 6*tau;
        double dh11 = 3*tau2 - 2*tau;
        double hk = t[k+1] - t[k];
        double dydt = (y[k]*dh00 + hk*d[k]*dh10 + y[k+1]*dh01 + hk*d[k+1]*dh11);
        return dydt / hk;  // chain rule: dτ/dt = 1/hk
    }

    // Second derivative d²y/dt²
    double deriv2(double s) const
    {
        double tau; int k = locate(s, tau);
        double ddh00 = 12*tau - 6;
        double ddh10 = 6*tau - 4;
        double ddh01 = -12*tau + 6;
        double ddh11 = 6*tau - 2;
        double hk = t[k+1] - t[k];
        double d2ydtau2 = y[k]*ddh00 + hk*d[k]*ddh10 + y[k+1]*ddh01 + hk*d[k+1]*ddh11;
        return d2ydtau2 / (hk * hk);
    }
};

// ── Degree-4 Lagrange interpolant for phi(t) at exactly 5 nodes ──────────
//
// phi(t) = sum_{i=0}^{4} phi[i] * L_i(t)  — unique degree-4 polynomial.
// C^∞ on ℝ: no interior knots, so phi''' is continuous everywhere.
// Compare: PCHIP is C^1 at knots, natural cubic spline is C^2 at knots.
//
// Uses the numerically stable barycentric second form:
//   phi(t) = [sum_i w_i*phi_i/(t-t_i)] / [sum_i w_i/(t-t_i)]
//   phi'(t) = [sum_i w_i*(phi(t)-phi_i)/(t-t_i)^2] / [sum_i w_i/(t-t_i)]
// At node t_k an exact limit formula avoids the 0/0 cancellation.
//
// Monotonicity is NOT guaranteed (unlike PCHIP).  Call is_monotone()
// before use and fall back to Pchip5 if it returns false.

struct LagrangePhi5 {
    double t_[5], y_[5], w_[5];  // nodes, values, barycentric weights

    void init(double const ts[5], double const ys[5])
    {
        for (int i = 0; i < 5; ++i) { t_[i] = ts[i]; y_[i] = ys[i]; }
        // w_i = 1 / prod_{j!=i} (t_i - t_j)
        for (int i = 0; i < 5; ++i) {
            w_[i] = 1.0;
            for (int j = 0; j < 5; ++j)
                if (j != i) w_[i] /= (ts[i] - ts[j]);
        }
    }

    // Barycentric second form — exact at nodes (guarded against 0/0).
    double eval(double s) const
    {
        double num = 0.0, den = 0.0;
        for (int i = 0; i < 5; ++i) {
            double d = s - t_[i];
            if (std::abs(d) < 1e-14 * (t_[4] - t_[0])) return y_[i];
            double tmp = w_[i] / d;
            num += tmp * y_[i];
            den += tmp;
        }
        return num / den;
    }

    // First derivative via barycentric differentiation.
    // Non-node: phi'(s) = [sum_i w_i*(phi(s)-y_i)/(s-t_i)^2] / [sum_i w_i/(s-t_i)]
    // At node t_k: phi'(t_k) = (1/w_k) * sum_{j!=k} w_j*(y_j-y_k)/(t_k-t_j)
    double deriv1(double s) const
    {
        double range = t_[4] - t_[0];
        for (int i = 0; i < 5; ++i) {
            if (std::abs(s - t_[i]) < 1e-14 * range) {
                double sum = 0.0;
                for (int j = 0; j < 5; ++j)
                    if (j != i) sum += w_[j] * (y_[j] - y_[i]) / (t_[i] - t_[j]);
                return sum / w_[i];
            }
        }
        double phi_s = eval(s);
        double num = 0.0, den = 0.0;
        for (int i = 0; i < 5; ++i) {
            double d = s - t_[i];
            num += w_[i] * (phi_s - y_[i]) / (d * d);
            den += w_[i] / d;
        }
        return num / den;
    }

    // Sample phi'(t) at N_SAMP interior points to verify strict monotonicity.
    // increasing=true → require phi'>0, false → phi'<0.
    bool is_monotone(bool increasing) const
    {
        constexpr int N_SAMP = 40;
        for (int k = 0; k <= N_SAMP; ++k) {
            double s = t_[0] + (t_[4] - t_[0]) * k / N_SAMP;
            double d = deriv1(s);
            if (increasing  && d <= 0.0) return false;
            if (!increasing && d >= 0.0) return false;
        }
        return true;
    }
};

// ── Best-fit 2D plane through n points in ℝ^Dim ─────────────────────────
//
// Returns: center (mean), e1 (direction of most variance), e2 (second).
// e1 is additionally oriented to point along (pts[n-1] - pts[0]) — forward.
// e2 = Gram-Schmidt of the second principal component w.r.t. e1.
//
// Algorithm: covariance matrix (Dim×Dim) → SymEig<Dim> → top 2 eigenvectors.

template <int Dim>
struct PlaneND {
    VecN<Dim> center;
    VecN<Dim> e1, e2;  // orthonormal in-plane basis
    double pts2d[5][2]; // projections of 5 input points
};

template <int Dim>
inline PlaneND<Dim> best_fit_plane(VecN<Dim> const pts[5])
{
    PlaneND<Dim> pl;

    // Mean
    pl.center = VecN<Dim>{};
    for (int i = 0; i < 5; ++i) pl.center += pts[i];
    pl.center = pl.center * (1.0/5.0);

    // Build Dim×Dim covariance matrix (as a flat array for SymEig)
    double cov[Dim][Dim] = {};
    for (int i = 0; i < 5; ++i) {
        VecN<Dim> q = pts[i] - pl.center;
        for (int r = 0; r < Dim; ++r)
            for (int c = 0; c < Dim; ++c)
                cov[r][c] += q[r] * q[c];
    }

    // Eigendecomposition of covariance
    SymEig<Dim> eig;
    eig.compute(cov);

    // Two largest eigenvectors (last two in ascending order)
    VecN<Dim> raw_e1{}, raw_e2{};
    for (int i = 0; i < Dim; ++i) {
        raw_e1[i] = eig.vec[Dim-1][i];  // largest eigenvalue → most variance
        raw_e2[i] = eig.vec[Dim-2][i];  // second largest
    }

    // Orient e1 along forward traversal direction
    VecN<Dim> chord = pts[4] - pts[0];
    if (chord.dot(raw_e1) < 0.0) raw_e1 = raw_e1 * -1.0;

    // Gram-Schmidt e2 perpendicular to e1 (raw_e2 may not be exactly orthogonal
    // after orientation flip, but SymEig produces orthogonal eigenvectors)
    raw_e2 = raw_e2 - raw_e2.dot(raw_e1) * raw_e1;
    double n2 = raw_e2.norm();
    if (n2 < 1e-15) {
        // Degenerate plane (all 5 points collinear) — pick any perpendicular
        // Use a simple construction: find component of e1 with smallest magnitude,
        // set that to 1, zero others, then Gram-Schmidt
        int mn = 0;
        for (int i = 1; i < Dim; ++i)
            if (std::abs(raw_e1[i]) < std::abs(raw_e1[mn])) mn = i;
        raw_e2 = VecN<Dim>{};
        raw_e2[mn] = 1.0;
        raw_e2 = raw_e2 - raw_e2.dot(raw_e1) * raw_e1;
        n2 = raw_e2.norm();
    }
    raw_e2 = raw_e2 * (1.0/n2);

    pl.e1 = raw_e1;
    pl.e2 = raw_e2;

    // Project 5 points onto plane
    for (int i = 0; i < 5; ++i) {
        VecN<Dim> q = pts[i] - pl.center;
        pl.pts2d[i][0] = q.dot(pl.e1);
        pl.pts2d[i][1] = q.dot(pl.e2);
    }

    return pl;
}

}  // namespace detail


// ── ConicWindow<Dim>: 5-point conic arc window ────────────────────────────
//
// Constructor: fits a conic to 5 control points (in their best-fit 2D plane),
// then builds a cross-ratio Lagrange interpolant (C^∞, affinely invariant).
// valid() is true iff the orbit is usable (r_i monotone, ctrl_err ≤ 1e-3).
//
// operator()(t) evaluates the conic arc at parameter t.
//
// If ConicWindow construction were exposed as a public API, the natural C++23
// form would be:
//
//   #if __cplusplus >= 202302L
//   template <int Dim>
//   std::expected<ConicWindow<Dim>, std::string>
//   make_conic_window(VecN<Dim> p0..p4, double t0..t4);
//   #endif
//
// That gives callers per-window failure info (which window failed and why)
// without exceptions, and composes cleanly with expected-based pipelines.
// The current blend_curve(conic_tag) all-or-nothing fallback is sufficient
// for the high-level API; the per-window form matters when callers want to
// inspect or override individual window choices.

template <int Dim>
class ConicWindow {
    bool valid_;

    // Plane parameters (stored for orbit evaluation)
    VecN<Dim> center_;  // 3D/nD centroid of the 5 control points
    VecN<Dim> e1_, e2_; // in-plane orthonormal basis

    // Principal frame rotation (2×2): principal-frame → SVD local frame
    // Stored as: V_pf[col][row], so point_local = pts2d @ V_pf
    double V_pf_[2][2];
    bool swapped_;    // true if x↔y coordinate swap was applied in principal frame

    // Rational map parameters (in principal frame, B=0)
    // P₀ = algebraic vertex (where L_e=0) when computable; inner-3 argmin |L/M| otherwise.
    double x0_, y0_;  // P₀ coordinates
    double A_e_, C_e_, L_e_, M_e_;  // L_e=2*A_e*x0+D_e (=0 at vertex), M_e=2*C_e*y0+E_e; B_e=0

    // phi(t) interpolant — stores φ_i = 2·arctan(s_i), unwrapped.
    // Lagrange quartic (C^∞) when monotone, else PCHIP (C^1).
    detail::LagrangePhi5 phi_lagrange_;
    detail::Pchip5       phi_pchip_;
    bool use_lagrange_;
    double ts_[5];   // control times (for range queries)

    // Line mode: collinear input is a valid degenerate conic (a straight line).
    // When active, eval_at_ uses degree-4 Lagrange in the plane instead of the
    // rational conic map.
    bool   line_mode_;
    double line_x_[5], line_y_[5];

    // Orbit reconstruction quality: max over the 5 control points of
    //   |eval(ts[i]) − projected_pt[i]| / win_scale
    // where win_scale = max pairwise distance of the 2D-projected points.
    // Dimensionless and invariant under rigid motions and uniform scaling.
    // 0.0 in line mode (degree-4 Lagrange interpolates exactly).
    double fit_error_;

public:
    ConicWindow() : valid_(false), line_mode_(false), fit_error_(0.0) {}

    // allow_cross_branch: if true, windows that span two hyperbola branches are
    // accepted using φ = 2·arctan(s) with unwrapping (projective arc through ∞).
    // Output may contain large/NaN values where the asymptote crosses the blend
    // region; the caller is responsible for handling or avoiding that range.
    ConicWindow(
        VecN<Dim> const& p0, VecN<Dim> const& p1, VecN<Dim> const& p2,
        VecN<Dim> const& p3, VecN<Dim> const& p4,
        double t0, double t1, double t2, double t3, double t4,
        bool allow_cross_branch = false)
        : valid_(false)
    {
        VecN<Dim> pts5[5] = {p0, p1, p2, p3, p4};
        ts_[0]=t0; ts_[1]=t1; ts_[2]=t2; ts_[3]=t3; ts_[4]=t4;

        // 1. Project to best-fit 2D plane
        detail::PlaneND<Dim> plane = detail::best_fit_plane<Dim>(pts5);
        center_ = plane.center;
        e1_     = plane.e1;
        e2_     = plane.e2;
        double (*pts2d)[2] = plane.pts2d;  // 5×2

        // 1b. Collinearity guard: if 5 projected 2D points are nearly collinear
        // (e.g., symmetric helix windows whose endpoints project to the same 2D point),
        // the conic fit is degenerate.  Check via 2D covariance: if the smaller
        // eigenvalue of the 2D scatter is < 1e-4 * larger, the points are nearly
        // collinear and the window is geometrically invalid for conic blending.
        {
            double mx=0,my=0;
            for(int i=0;i<5;++i){mx+=pts2d[i][0];my+=pts2d[i][1];}
            mx/=5; my/=5;
            double sxx=0,sxy=0,syy=0;
            for(int i=0;i<5;++i){
                double dx=pts2d[i][0]-mx, dy=pts2d[i][1]-my;
                sxx+=dx*dx; sxy+=dx*dy; syy+=dy*dy;
            }
            // eigenvalues of [[sxx,sxy],[sxy,syy]]
            double tr=sxx+syy, dsc=(sxx-syy)*(sxx-syy)+4*sxy*sxy;
            double lam_min = 0.5*(tr - std::sqrt(dsc));
            double lam_max = 0.5*(tr + std::sqrt(dsc));
            if (lam_max < 1e-30 || lam_min < 1e-4 * lam_max) {
                // Straight line: a valid degenerate conic — monotone by definition.
                // Store 2D projections; eval_at_ will use degree-4 Lagrange.
                for (int i = 0; i < 5; ++i) {
                    line_x_[i] = pts2d[i][0];
                    line_y_[i] = pts2d[i][1];
                }
                line_mode_  = true;
                fit_error_  = 0.0;   // degree-4 Lagrange interpolates exactly
                valid_      = true;
                return;
            }
        }

        // 2. Fit conic coefficients [A,B,C,D,E,F]
        double coeffs[6];
        detail::fit_conic_2d(pts2d, coeffs);
        double A = coeffs[0], B = coeffs[1], C = coeffs[2];
        double D = coeffs[3], E = coeffs[4], F = coeffs[5];
        (void)F;

        // 3. Rotate to principal frame (B_p=0 by construction)
        detail::Eig2x2 epf = detail::eigh2(A, B*0.5, C);

        // Sort by |eigenvalue| ascending → consistent assignment
        int i0 = (std::abs(epf.val[0]) <= std::abs(epf.val[1])) ? 0 : 1;
        int i1 = 1 - i0;
        double lam0 = epf.val[i0], lam1 = epf.val[i1];
        double vx0 = epf.vec[i0][0], vy0 = epf.vec[i0][1];
        double vx1 = epf.vec[i1][0], vy1 = epf.vec[i1][1];

        // Sign convention: first principal axis has positive global x-component
        // (dot with e1 in the original space)
        double phys_x = vx0 * e1_[0] + vy0 * e2_[0];
        double phys_y = vx0 * e1_[1] + vy0 * e2_[1];
        if (phys_x < 0.0 || (std::abs(phys_x) < 1e-12 && phys_y < 0.0)) {
            vx0 = -vx0;  vy0 = -vy0;
        }
        // Maintain right-hand orientation (det > 0)
        if (vx0*vy1 - vy0*vx1 < 0.0) {
            vx1 = -vx1;  vy1 = -vy1;
        }

        // V_pf: columns are principal-frame axes in SVD-local coordinates
        // pts_p[i] = [pts2d[i][0]*vx0 + pts2d[i][1]*vy0,
        //             pts2d[i][0]*vx1 + pts2d[i][1]*vy1]
        V_pf_[0][0] = vx0;  V_pf_[0][1] = vx1;
        V_pf_[1][0] = vy0;  V_pf_[1][1] = vy1;

        // Project to principal frame
        double pts_p[5][2];
        for (int i = 0; i < 5; ++i) {
            pts_p[i][0] = pts2d[i][0]*vx0 + pts2d[i][1]*vy0;
            pts_p[i][1] = pts2d[i][0]*vx1 + pts2d[i][1]*vy1;
        }

        // Principal-frame conic coefficients (B_p=0)
        double A_p = lam0, C_p = lam1, B_p = 0.0;
        double DE[2] = {D, E};
        double D_p = DE[0]*vx0 + DE[1]*vy0;
        double E_p = DE[0]*vx1 + DE[1]*vy1;

        // 4. Swap x↔y if A_p near zero (near-parabola with vertical axis)
        swapped_ = false;
        double A_e = A_p, C_e = C_p, D_e = D_p, E_e = E_p;
        double pts_e[5][2];
        for (int i = 0; i < 5; ++i) {
            pts_e[i][0] = pts_p[i][0];
            pts_e[i][1] = pts_p[i][1];
        }

        constexpr double EPS_SWAP = 1e-9;
        if (std::abs(A_e) < EPS_SWAP * std::max(std::abs(C_e), 1.0)) {
            // Swap x↔y
            std::swap(A_e, C_e);
            std::swap(D_e, E_e);
            for (int i = 0; i < 5; ++i) std::swap(pts_e[i][0], pts_e[i][1]);
            swapped_ = true;
        }

        // 4b. Cross-branch detection (needed when allow_cross_branch is set).
        //
        // A window is "cross-branch" when consecutive control points lie on opposite
        // branches of a hyperbola.  Detection: project each pt (relative to conic
        // center) onto the transverse eigenvector; sign changes indicate branch crossings.
        //
        // lam0, lam1 are already the principal-frame eigenvalues (smaller |val| first).
        // is_hyperbola ↔ they have opposite signs ↔ 4AC − B² < 0.
        bool is_cross_branch = false;
        // Near-parabola guard: if |lam0| << |lam1| the smaller eigenvalue is
        // numerically zero (true parabola).  A tiny negative lam0 would give
        // lam0*lam1 < 0 and trigger the cross-branch check, but the "center"
        // then lies at infinity → catastrophic cancellation in the sign test.
        // Skip cross-branch detection for near-parabola windows.
        if (lam0 * lam1 < 0.0 &&
            std::abs(lam0) > 1e-6 * std::abs(lam1)) {
            // Conic center in SVD local frame
            double det4 = 4.0*A*C - B*B;  // < 0 for hyperbola
            double cx = (-2.0*C*D + B*E) / det4;
            double cy = (-2.0*A*E + B*D) / det4;
            // Transverse eigenvector = eigenvector whose eigenvalue has the
            // OPPOSITE sign to F.  For A_p·U² + C_p·V² + F = 0 the real
            // points satisfy A_p·U² = −F, so the transverse eigenvalue is
            // the one with sign(lam) ≠ sign(F).  F is invariant under the
            // rotation to principal frame and fit_conic_2d returns an
            // arbitrary-sign null vector, so we cannot rely on "positive
            // eigenvalue = transverse axis".
            double tx = (lam0 * F < 0.0) ? vx0 : vx1;
            double ty = (lam0 * F < 0.0) ? vy0 : vy1;
            for (int i = 0; i < 4 && !is_cross_branch; ++i) {
                double s0p = (pts2d[i][0]   - cx)*tx + (pts2d[i][1]   - cy)*ty;
                double s1p = (pts2d[i+1][0] - cx)*tx + (pts2d[i+1][1] - cy)*ty;
                if (s0p * s1p < 0.0) is_cross_branch = true;
            }
        }

        // 5–6. Find P₀ for φ = 2·arctan(s) parametrization (mirrors conicspline.py).
        //
        // Strategy:
        //   Primary: algebraic vertex (where L_e = 2*A_e*x + D_e = 0, so s₀=0).
        //     This is an intrinsic conic property — same vertex for any window on
        //     the same conic → identical orbits for both overlapping windows →
        //     exact blend for locally-conic data.  For parabolas (C_e≈0) the vertex
        //     exists and gives trivially monotone φ (s_i = x_i, arctan is monotone).
        //   Fallback: argmin |L_k/M_k| over inner-3 control points (k=1,2,3).
        //     These three points are shared by both overlapping windows at every
        //     blend junction, so the same argmin is chosen by both → same P₀.
        //
        // Parametrization: φ_i = 2·arctan(s_i), s_i = (y_i−y₀)/(x_i−x₀).
        // After ±π unwrapping, φ is monotone iff the arc doesn't double back.
        // Works uniformly for same-branch and cross-branch windows; no separate
        // cross-ratio path needed.

        constexpr double EPS_COEFF  = 1e-12;
        constexpr double EPS_DX     = 1e-9;
        constexpr double EPS_VERTEX = 1e-9;
        constexpr double PI         = 3.14159265358979323846;

        // Pre-compute L_k, M_k (conic gradient components) at all 5 control points
        double L_all[5], M_all[5];
        for (int k = 0; k < 5; ++k) {
            L_all[k] = 2.0*A_e*pts_e[k][0] + D_e;
            M_all[k] = 2.0*C_e*pts_e[k][1] + E_e;
        }

        // --- Step 5a: Try algebraic vertex as P₀ ---
        // Solve L_e = 0: x_v = -D_e/(2*A_e).
        // Then solve the conic for y_v; pick root closest to control-point centroid.
        bool   has_vertex = false;
        double x0_vert = 0.0, y0_vert = 0.0;

        if (std::abs(A_e) >= EPS_VERTEX * std::max(std::abs(C_e), 1.0)) {
            double x_v      = -D_e / (2.0 * A_e);
            double cv_const = A_e*x_v*x_v + D_e*x_v + F;  // = F − D_e²/(4*A_e)
            if (std::abs(C_e) < EPS_VERTEX) {
                // Parabola-like: linear in y → E_e*y + cv_const = 0
                if (std::abs(E_e) >= EPS_VERTEX) {
                    double y_v = -cv_const / E_e;
                    if (std::abs(2.0*C_e*y_v + E_e) >= EPS_COEFF) {
                        x0_vert = x_v; y0_vert = y_v; has_vertex = true;
                    }
                }
            } else {
                double disc = E_e*E_e - 4.0*C_e*cv_const;
                if (disc >= 0.0) {
                    double sq = std::sqrt(disc);
                    double y1 = (-E_e + sq) / (2.0*C_e);
                    double y2 = (-E_e - sq) / (2.0*C_e);
                    double cy_mean = 0.0;
                    for (int i = 0; i < 5; ++i) cy_mean += pts_e[i][1];
                    cy_mean /= 5.0;
                    double y_v = (std::abs(y1-cy_mean) <= std::abs(y2-cy_mean)) ? y1 : y2;
                    if (std::abs(2.0*C_e*y_v + E_e) >= EPS_COEFF) {
                        x0_vert = x_v; y0_vert = y_v; has_vertex = true;
                    }
                }
            }
        }

        // --- Step 5b: Helper to compute φ_i = 2·arctan(s_i) and check monotonicity ---
        // k_self: index of P₀ in pts_e (use tangent slope there), or -1 if vertex.
        auto phi_mono = [&](double x0, double y0, double Le, double Me,
                            int k_self, double pv[5]) -> bool
        {
            double sv[5];
            for (int i = 0; i < 5; ++i) {
                if (i == k_self) {
                    sv[i] = -Le / Me;   // tangent slope at control-point P₀
                } else {
                    double dx = pts_e[i][0] - x0;
                    double dy = pts_e[i][1] - y0;
                    if (std::abs(dx) >= EPS_DX * (std::abs(dy) + 1.0))
                        sv[i] = dy / dx;
                    else
                        sv[i] = (dy >= 0.0 ? 1.0 : -1.0) * 1e15;
                }
            }
            for (int i = 0; i < 5; ++i)
                pv[i] = 2.0 * std::atan(sv[i]);
            for (int i = 0; i < 4; ++i) {
                double d = pv[i+1] - pv[i];
                while (d >  PI) { d -= 2.0*PI; pv[i+1] -= 2.0*PI; }
                while (d < -PI) { d += 2.0*PI; pv[i+1] += 2.0*PI; }
            }
            bool mono_up = true, mono_dn = true;
            for (int i = 0; i < 4; ++i) {
                if (pv[i+1] <= pv[i]) mono_up = false;
                if (pv[i+1] >= pv[i]) mono_dn = false;
            }
            return mono_up || mono_dn;
        };

        // --- Step 5c: Try P₀ candidates, accept first with monotone φ ---
        bool   found    = false;
        double phi_vals[5];
        double x0_used  = 0.0, y0_used  = 0.0;
        double L_e_used = 0.0, M_e_used = 0.0;

        if (has_vertex) {
            double M_v = 2.0*C_e*y0_vert + E_e;
            if (phi_mono(x0_vert, y0_vert, 0.0, M_v, /*k_self=*/-1, phi_vals)) {
                x0_used = x0_vert; y0_used = y0_vert;
                L_e_used = 0.0;    M_e_used = M_v;
                found = true;
            }
        }

        if (!found) {
            // Argmin |L/M| over inner-3 (indices 1,2,3)
            double best = 1e30;
            int k_best = -1;
            for (int k = 1; k <= 3; ++k) {
                double Mk = M_all[k];
                if (std::abs(Mk) < EPS_COEFF) continue;
                double r = std::abs(L_all[k] / Mk);
                if (r < best) { best = r; k_best = k; }
            }
            if (k_best >= 0) {
                double Lk = L_all[k_best], Mk = M_all[k_best];
                if (phi_mono(pts_e[k_best][0], pts_e[k_best][1], Lk, Mk, k_best, phi_vals)) {
                    x0_used = pts_e[k_best][0]; y0_used = pts_e[k_best][1];
                    L_e_used = Lk; M_e_used = Mk;
                    found = true;
                }
            }
        }

        // Cross-branch windows are only accepted when caller explicitly permits it.
        if (found && is_cross_branch && !allow_cross_branch) found = false;

        if (!found) return;  // no P₀ gave monotone φ → invalid

        x0_ = x0_used; y0_ = y0_used;
        A_e_ = A_e;  C_e_ = C_e;  L_e_ = L_e_used;  M_e_ = M_e_used;

        // 7. Build phi(t) interpolant.
        //    Try the degree-4 Lagrange polynomial first: it is C^∞ everywhere
        //    (no interior knots), so phi''' is continuous at ts[2] = the blend
        //    junction, giving C^3 of the orbit there → τ continuity.
        //    Fall back to PCHIP (C^1) only if the Lagrange polynomial is
        //    non-monotone (which can happen for rapidly-varying phi data).
        {
            bool inc = (phi_vals[4] > phi_vals[0]);
            phi_lagrange_.init(ts_, phi_vals);
            if (phi_lagrange_.is_monotone(inc)) {
                use_lagrange_ = true;
            } else {
                phi_pchip_.init(ts_, phi_vals);
                use_lagrange_ = false;
            }
        }

        // 8. Validate: check orbit reconstructs all 5 projected control points.
        // eval_at_ always returns a point in the best-fit plane, so compare against
        // the plane projections (center_ + pts2d[i][0]*e1_ + pts2d[i][1]*e2_), NOT
        // the original nD pts5[i].  For a non-planar curve (helix, torus, …) pts5[i]
        // has a planarity residual ⊥ to the plane that the orbit can never reproduce.
        //
        // Normalise by win_scale = max pairwise distance of the 2D-projected points
        // so the threshold is dimensionless and invariant under rigid motions + scaling.
        double wsq = 0.0;
        for (int i = 0; i < 5; ++i)
            for (int j = i+1; j < 5; ++j) {
                double dx = pts2d[i][0]-pts2d[j][0], dy = pts2d[i][1]-pts2d[j][1];
                double d2 = dx*dx + dy*dy;
                if (d2 > wsq) wsq = d2;
            }
        double win_scale = std::sqrt(wsq);

        double max_err = 0.0;
        for (int i = 0; i < 5; ++i) {
            VecN<Dim> pred     = eval_at_(ts_[i]);
            VecN<Dim> proj_pt  = center_ + pts2d[i][0] * e1_ + pts2d[i][1] * e2_;
            double err = (pred - proj_pt).norm();
            if (err > max_err) max_err = err;
        }
        fit_error_ = (win_scale > 0.0) ? max_err / win_scale : 0.0;
        if (fit_error_ > 1e-3) return;  // conic orbit doesn't reconstruct projected control points

        valid_ = true;
    }

    bool   valid()        const { return valid_; }
    double fit_error()    const { return fit_error_; }   // dimensionless; 0 = exact
    bool   uses_c_infty() const { return use_lagrange_; }  // true → Lagrange (C^∞ phi)

    VecN<Dim> operator()(double t) const { return eval_at_(t); }

private:
    // Rational map: u(s) = -(L_e + M_e*s)/Q(s),  x=x0+u,  y=y0+s*u
    // In principal frame (B_e=0):
    //   Q(s) = A_e + C_e*s²
    //   L_e = 2*A_e*x0 + D_e  (gradient x at P0; ≈0 for canonical P0)
    //   u(s) = -(L_e + M_e*s) / (A_e + C_e*s²)
    VecN<Dim> eval_at_(double t) const
    {
        // Line mode: degenerate conic (collinear input).  Degree-4 Lagrange in 2D.
        if (line_mode_) {
            double x = 0.0, y = 0.0;
            for (int j = 0; j < 5; ++j) {
                double L = 1.0;
                for (int k = 0; k < 5; ++k)
                    if (k != j) L *= (t - ts_[k]) / (ts_[j] - ts_[k]);
                x += L * line_x_[j];
                y += L * line_y_[j];
            }
            return center_ + x * e1_ + y * e2_;
        }

        // φ = 2·arctan(s) → s = tan(φ/2).  Near φ = ±π the arc crosses ∞ (asymptote).
        double param = use_lagrange_ ? phi_lagrange_.eval(t) : phi_pchip_.eval(t);
        double s = std::tan(param * 0.5);

        double Q = A_e_ + C_e_ * s * s;
        double x2d, y2d;
        constexpr double EPS_DET = 1e-12;
        if (std::abs(Q) < EPS_DET) {
            // Asymptote: return NaN-like (large value won't pass ctrl_err check)
            x2d = y2d = 1e300;
        } else {
            double u = -(L_e_ + M_e_ * s) / Q;
            x2d = x0_ + u;
            y2d = y0_ + s * u;
        }

        // Swap back if needed
        if (swapped_) std::swap(x2d, y2d);

        // Principal frame → SVD local frame: [lx, ly] = [x2d, y2d] @ V_pf^T
        // V_pf stores M = [[vx0,vx1],[vy0,vy1]] (row-major, rows = SVD coords).
        // Forward: [x_p,y_p] = [x_svd,y_svd] @ M
        // Inverse: [lx,ly]   = [x_p,y_p]     @ M^T  (M orthogonal → M^{-1}=M^T)
        // M^T[0,0]=vx0=V_pf_[0][0], M^T[1,0]=vx1=V_pf_[0][1]
        // M^T[0,1]=vy0=V_pf_[1][0], M^T[1,1]=vy1=V_pf_[1][1]
        double lx = x2d * V_pf_[0][0] + y2d * V_pf_[0][1];
        double ly = x2d * V_pf_[1][0] + y2d * V_pf_[1][1];

        // SVD local frame → nD space
        return center_ + lx * e1_ + ly * e2_;
    }
};

// ── LagrangeWindow: per-window fallback for invalid ConicWindows ────────────
//
// Degree-4 Lagrange polynomial interpolant through 5 control points.
// Used in place of a ConicWindow when the conic fit fails (degenerate
// geometry, non-monotone cross-ratio, etc.).
//
// Always valid(); evaluation is O(25) multiplications per call.
// C^∞ within its parameter range; C^N at junctions because each junction
// is covered by exactly one window used identically from both adjacent segments.

template <int Dim>
struct LagrangeWindow {
    VecN<Dim> pts_[5];
    double    ts_[5];

    LagrangeWindow() = default;

    LagrangeWindow(VecN<Dim> const& p0, VecN<Dim> const& p1,
                   VecN<Dim> const& p2, VecN<Dim> const& p3,
                   VecN<Dim> const& p4,
                   double t0, double t1, double t2, double t3, double t4)
    {
        pts_[0]=p0; pts_[1]=p1; pts_[2]=p2; pts_[3]=p3; pts_[4]=p4;
        ts_[0]=t0;  ts_[1]=t1;  ts_[2]=t2;  ts_[3]=t3;  ts_[4]=t4;
    }

    bool valid()        const { return true; }
    bool uses_c_infty() const { return true; }

    VecN<Dim> operator()(double t) const {
        VecN<Dim> r{};
        for (int i = 0; i < 5; ++i) {
            double L = 1.0;
            for (int j = 0; j < 5; ++j)
                if (j != i) L *= (t - ts_[j]) / (ts_[i] - ts_[j]);
            r = r + pts_[i] * L;
        }
        return r;
    }
};

// ── Tagged overload: blend_curve(..., conic_tag{}) ─────────────────────────
//
// Builds ConicWindow objects for each group of 5 consecutive control points.
// Valid windows use the conic cross-ratio parametrisation; invalid windows
// fall back to a LagrangeWindow (degree-4 polynomial) for that window only.
// No global fallback: a local invalidity has only local effect.
//
// Minimum control points: 6  (gives 1 blended segment j=2).
// Blended segments: j = 2, …, n-4.
// At junction times[j]: both adjacent segments use win[j-2] → C^N continuity.
//
// used_conic (optional out-parameter): true if at least one conic window was
// used, false if all windows fell back to Lagrange.  Pass nullptr to ignore.

template <int Dim>
inline BlendResultND<Dim> blend_curve(
#if __cplusplus >= 202002L
    std::span<VecN<Dim> const> ctrl,
    std::span<double const>    times,
#else
    std::vector<VecN<Dim>> const& ctrl,
    std::vector<double>    const& times,
#endif
    conic_tag,
    int   pts_per_seg        = 60,
    int   smooth_N           = 2,
    bool* used_conic         = nullptr,
    bool  allow_cross_branch = false)
{
    int n = static_cast<int>(ctrl.size());
    if (n != static_cast<int>(times.size()))
        FC_THROW("fc::blend_curve(conic): ctrl and times must have equal length");
    if (n < 6)
        FC_THROW("fc::blend_curve(conic): need at least 6 control points");
    if (pts_per_seg < 2)
        FC_THROW("fc::blend_curve(conic): pts_per_seg must be >= 2");

    // Build (n-4) windows: ConicWindow where valid, LagrangeWindow otherwise.
    // Per-window fallback: invalidity is local — no global effect on other windows.
    using Win = std::variant<ConicWindow<Dim>, LagrangeWindow<Dim>>;
    std::vector<Win> wins;
    wins.reserve(n - 4);
    int n_conic = 0;

    for (int i = 0; i < n-4; ++i) {
        ConicWindow<Dim> cw(ctrl[i], ctrl[i+1], ctrl[i+2], ctrl[i+3], ctrl[i+4],
                            times[i], times[i+1], times[i+2], times[i+3], times[i+4],
                            allow_cross_branch);
        if (cw.valid()) {
            ++n_conic;
            wins.emplace_back(std::move(cw));
        } else {
            wins.emplace_back(LagrangeWindow<Dim>(
                ctrl[i], ctrl[i+1], ctrl[i+2], ctrl[i+3], ctrl[i+4],
                times[i], times[i+1], times[i+2], times[i+3], times[i+4]));
        }
    }
    if (used_conic) *used_conic = (n_conic > 0);

    // Evaluate any window type through std::visit.
    auto eval_win = [](Win const& w, double t) -> VecN<Dim> {
        return std::visit([t](auto const& ww) -> VecN<Dim> { return ww(t); }, w);
    };

    // Blend: segment j runs from times[j] to times[j+1], j=2..n-4
    // Uses win[j-2] (w=0 side) and win[j-1] (w=1 side)
    // At s=0: blend = win[j-2](times[j])  = ctrl[j]   ✓
    // At s=1: blend = win[j-1](times[j+1]) = ctrl[j+1] ✓
    int n_segs = n - 5;
    BlendResultND<Dim> out;
    out.pts.reserve(n_segs * pts_per_seg + 1);
    out.times.reserve(n_segs * pts_per_seg + 1);

    for (int j = 2; j <= n-4; ++j) {
        Win const& wA = wins[j-2];
        Win const& wB = wins[j-1];
        double t_lo = times[j], t_hi = times[j+1];

        int k_max = (j < n-4) ? pts_per_seg - 1 : pts_per_seg;
        for (int k = 0; k <= k_max; ++k) {
            double s = static_cast<double>(k) / pts_per_seg;
            double t = t_lo + s * (t_hi - t_lo);
            double w = smoothstep(s, smooth_N);
            VecN<Dim> A = eval_win(wA, t), B = eval_win(wB, t);
            out.pts.push_back(A*(1.0-w) + B*w);
            out.times.push_back(t);
        }
    }

    return out;
}

// ── 3D non-template overloads ─────────────────────────────────────────────

inline BlendResult blend_curve(
#if __cplusplus >= 202002L
    std::span<Vec3 const>   ctrl,
    std::span<double const> times,
#else
    std::vector<Vec3>   const& ctrl,
    std::vector<double> const& times,
#endif
    conic_tag,
    int   pts_per_seg        = 60,
    int   smooth_N           = 2,
    bool* used_conic         = nullptr,
    bool  allow_cross_branch = false)
{
    return blend_curve<3>(ctrl, times, conic_tag{}, pts_per_seg, smooth_N, used_conic, allow_cross_branch);
}

}  // namespace fc
