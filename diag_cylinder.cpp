// diag_cylinder.cpp — does an exact right circular cylinder exist through 5 twisted cubic pts?
//
// For each 5-point window on the twisted cubic r(t)=(t,t²,t³):
//   1. Frenet shortcut: use e3=e1×e2 as fixed axis → circle-fit residual (NOT zero → not exact)
//   2. Optimized axis: minimize circle-fit residual over û ∈ S² (2D gradient descent)
//      If optimized residual → machine epsilon, an exact cylinder exists.
//      If optimized residual >> machine epsilon, no exact right circular cylinder passes
//      through those 5 points.
//
// The optimizer uses a 40×20 grid search + gradient descent refinement.
//
// Compile:
//   g++ -std=c++17 -O2 -I. -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 -o diag_cylinder diag_cylinder.cpp

#include "conicblend.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace fc;
using namespace fc::detail;
using V3 = VecN<3>;

// ── Twisted cubic ──────────────────────────────────────────────────────────
static V3 tc(double t) { return V3(t, t*t, t*t*t); }

// ── Cross product ──────────────────────────────────────────────────────────
static V3 cross3(V3 const& a, V3 const& b) {
    return V3(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}

// ── 2D circle fit (Bookstein algebraic, least-squares) ─────────────────────
// Fits x² + y² + Dx + Ey + F = 0.  cx = -D/2, cy = -E/2, r = sqrt(cx²+cy²-F).
// Works even when points are overdetermined (5 pts, 3 params → least squares).
static void fit_circle_2d(double pts[][2], int n,
                          double& cx, double& cy, double& r)
{
    double sx=0,sy=0,sxx=0,sxy=0,syy=0,sb=0,sbx=0,sby=0;
    for (int i=0;i<n;++i) {
        double x=pts[i][0], y=pts[i][1], b=x*x+y*y;
        sx+=x; sy+=y; sxx+=x*x; sxy+=x*y; syy+=y*y;
        sb+=b; sbx+=b*x; sby+=b*y;
    }
    double N=n;
    double A3[3][3]={{sxx,sxy,sx},{sxy,syy,sy},{sx,sy,N}};
    double b3[3]={-sbx,-sby,-sb};  // circle: Dx+Ey+F = -(x²+y²)
    for (int col=0;col<3;++col) {
        int pm=col;
        for (int row=col+1;row<3;++row)
            if (std::abs(A3[row][col])>std::abs(A3[pm][col])) pm=row;
        std::swap(A3[col],A3[pm]); std::swap(b3[col],b3[pm]);
        if (std::abs(A3[col][col])<1e-20) continue;
        double inv=1.0/A3[col][col];
        for (int row=col+1;row<3;++row) {
            double f=A3[row][col]*inv;
            for (int c=col;c<3;++c) A3[row][c]-=f*A3[col][c];
            b3[row]-=f*b3[col];
        }
    }
    double D_c=0,E_c=0,F_c=0;
    if (std::abs(A3[2][2])>1e-20) F_c=b3[2]/A3[2][2];
    if (std::abs(A3[1][1])>1e-20) E_c=(b3[1]-A3[1][2]*F_c)/A3[1][1];
    if (std::abs(A3[0][0])>1e-20) D_c=(b3[0]-A3[0][1]*E_c-A3[0][2]*F_c)/A3[0][0];
    cx=-0.5*D_c; cy=-0.5*E_c;
    double r2=cx*cx+cy*cy-F_c;
    r=(r2>=0)?std::sqrt(r2):0.0;
}

// ── RMS residual from circle fit for a given axis direction û ───────────────
static double cylinder_residual(V3 const pts5[5], V3 const u_hat)
{
    // Build orthonormal (v_hat, w_hat) perpendicular to u_hat
    V3 tmp = (std::abs(u_hat[0]) < 0.9) ? V3(1,0,0) : V3(0,1,0);
    V3 v_hat = tmp - tmp.dot(u_hat)*u_hat;
    double vn = v_hat.norm();
    if (vn < 1e-12) return 1e30;
    v_hat = v_hat * (1.0/vn);
    V3 w_hat = cross3(u_hat, v_hat);

    // Subtract centroid for numerical stability
    V3 cen{};
    for (int i=0;i<5;i++) cen = cen + pts5[i];
    cen = cen * 0.2;

    double disk[5][2];
    for (int i=0;i<5;i++) {
        V3 q = pts5[i] - cen;
        disk[i][0] = q.dot(v_hat);
        disk[i][1] = q.dot(w_hat);
    }

    double cx,cy,r;
    fit_circle_2d(disk, 5, cx, cy, r);

    double res=0.0;
    for (int i=0;i<5;i++) {
        double dx=disk[i][0]-cx, dy=disk[i][1]-cy;
        double d=std::sqrt(dx*dx+dy*dy)-r;
        res+=d*d;
    }
    return std::sqrt(res/5.0);
}

// ── Optimize û over S² to minimise cylinder_residual ──────────────────────
// Returns best residual; sets best_u.
static double optimize_cylinder(V3 const pts5[5], V3& best_u)
{
    // Phase 1: grid search on the sphere + explicit poles
    double best_res = 1e30;
    double best_th = M_PI/2, best_ph = 0;

    // Poles: (0,0,±1), (±1,0,0), (0,±1,0) — never in the sin/cos grid
    V3 poles[6]={V3(1,0,0),V3(-1,0,0),V3(0,1,0),V3(0,-1,0),V3(0,0,1),V3(0,0,-1)};
    for (auto& u : poles) {
        double r = cylinder_residual(pts5, u);
        if (r < best_res) {
            best_res=r;
            // back-convert to (th, ph)
            double th0=std::acos(std::max(-1.0,std::min(1.0,u[2])));
            double ph0=std::atan2(u[1],u[0]);
            best_th=th0; best_ph=ph0;
        }
    }

    int NTH=40, NPH=40;
    for (int it=0;it<NTH;++it) {
        double th = M_PI*(it+0.5)/NTH;
        for (int ip=0;ip<NPH;++ip) {
            double ph = 2*M_PI*ip/NPH;
            V3 u(std::sin(th)*std::cos(ph),
                 std::sin(th)*std::sin(ph),
                 std::cos(th));
            double r = cylinder_residual(pts5, u);
            if (r < best_res) { best_res=r; best_th=th; best_ph=ph; }
        }
    }

    // Phase 2: gradient descent with backtracking line search
    double th=best_th, ph=best_ph;
    double step=0.05;
    for (int iter=0; iter<2000 && step>1e-12; ++iter) {
        V3 u0(std::sin(th)*std::cos(ph), std::sin(th)*std::sin(ph), std::cos(th));
        double r0 = cylinder_residual(pts5, u0);

        // Finite-difference gradient in (th, ph)
        double eps=1e-6;
        double dth = (cylinder_residual(pts5,
                        V3(std::sin(th+eps)*std::cos(ph),
                           std::sin(th+eps)*std::sin(ph),
                           std::cos(th+eps))) - r0) / eps;
        double dph = (cylinder_residual(pts5,
                        V3(std::sin(th)*std::cos(ph+eps),
                           std::sin(th)*std::sin(ph+eps),
                           std::cos(th))) - r0) / eps;

        double gn = std::sqrt(dth*dth + dph*dph);
        if (gn < 1e-20) break;

        // Normalised gradient step
        double th2 = th - step*dth/gn;
        double ph2 = ph - step*dph/gn;
        V3 u2(std::sin(th2)*std::cos(ph2), std::sin(th2)*std::sin(ph2), std::cos(th2));
        double r2 = cylinder_residual(pts5, u2);

        if (r2 < r0) {
            th=th2; ph=ph2; best_res=r2;
            step *= 1.2;
        } else {
            step *= 0.5;
        }
    }

    best_u = V3(std::sin(th)*std::cos(ph), std::sin(th)*std::sin(ph), std::cos(th));
    return best_res;
}

// ── Analyse one 5-point window ──────────────────────────────────────────────
static void analyse_window(int win, V3 pts5[5])
{
    // ── Frenet shortcut residual ───────────────────────────────────────────
    PlaneND<3> plane = best_fit_plane<3>(pts5);
    V3 e1=plane.e1, e2=plane.e2;
    V3 e3=cross3(e1,e2);

    double plan_rms = 0.0;
    for (int i=0;i<5;++i) {
        V3 proj = plane.center + plane.pts2d[i][0]*e1 + plane.pts2d[i][1]*e2;
        double d2 = (pts5[i]-proj).dot(pts5[i]-proj);
        plan_rms += d2;
    }
    plan_rms = std::sqrt(plan_rms/5.0);

    double res_frenet = cylinder_residual(pts5, e3);

    // ── Optimized axis residual ────────────────────────────────────────────
    V3 best_u;
    double res_opt = optimize_cylinder(pts5, best_u);

    // Angle between e3 and optimized axis
    double cos_ang = std::abs(e3.dot(best_u));
    if (cos_ang > 1.0) cos_ang = 1.0;
    double ang_deg = std::acos(cos_ang) * 180.0 / M_PI;

    std::printf("win %-2d  planarity=%.2e  "
                "frenet_res=%.2e  opt_res=%.2e  axis_tilt=%.1f°  "
                "improvement=%5.1fx\n",
                win, plan_rms, res_frenet, res_opt, ang_deg,
                res_frenet / (res_opt + 1e-20));
}

// ── Helix: r(t) = (R cos ωt,  R sin ωt,  p·t) ──────────────────────────────
// Lies exactly on the right circular cylinder x²+y² = R², axis = (0,0,1).
static V3 helix(double t, double R, double omega, double pitch)
{
    return V3(R*std::cos(omega*t), R*std::sin(omega*t), pitch*t);
}

int main()
{
    std::printf("planarity   = RMS distance from best-fit plane (current approach's knot error)\n");
    std::printf("frenet_res  = cylinder residual with axis=e3 (Frenet shortcut)\n");
    std::printf("opt_res     = cylinder residual with optimized axis\n");
    std::printf("axis_tilt   = angle between e3 and optimized axis\n");
    std::printf("If opt_res → machine-eps (~1e-15), an exact cylinder exists.\n");
    std::printf("If opt_res >> machine-eps, no exact right circular cylinder exists.\n\n");
    std::printf("%-68s\n", std::string(68, '-').c_str());

    // ── Test A: helix (known exact cylinder) — validate the optimizer ─────
    {
        int n = 12;
        double R=1.0, omega=2.0, pitch=0.5;
        std::vector<V3> ctrl(n);
        std::vector<double> ts(n);
        for (int i=0;i<n;++i) {
            double t = 2.0*M_PI*(double)i/(n-1);
            ctrl[i] = helix(t, R, omega, pitch);
            ts[i] = t;
        }
        std::printf("\nTest A — HELIX  R=%.1f  ω=%.1f  pitch=%.1f  n=%d\n"
                    "  Exact cylinder axis = (0,0,1).  opt_res should → machine-eps.\n\n",
                    R, omega, pitch, n);
        for (int i=0;i<n-4;++i) {
            V3 pts5[5]={ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4]};
            analyse_window(i, pts5);
        }
    }

    // ── Test B: twisted cubic (no exact cylinder expected) ────────────────
    {
        int n = 12;
        std::vector<V3> ctrl(n);
        std::vector<double> ts(n);
        for (int i=0;i<n;++i) {
            double t=(double)i/(n-1);
            ctrl[i]=tc(t); ts[i]=t;
        }
        std::printf("\nTest B — TWISTED CUBIC  r(t)=(t,t²,t³)  n=%d\n"
                    "  opt_res > planarity → no exact cylinder; opt_res → planarity would be needed.\n\n",
                    n);
        for (int i=0;i<n-4;++i) {
            V3 pts5[5]={ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4]};
            analyse_window(i, pts5);
        }
    }
    return 0;
}
