// conicblend_cylinder.hpp — C^N 3D curve interpolation via cylinder unrolling
// Header-only C++17 library.  Requires conicblend.hpp.
//
// Architecture:
//   Each 5-point window finds the best-fitting right circular cylinder through
//   the 5 control points (smallest Δφ arc span), unrolls them to (r·φ, z), fits
//   a 2D ConicWindow there, and re-rolls the result back to 3D.
//
//   Fallback chain per window:
//     CylinderWindow (cylinder found, 2D conic valid)
//     → ConicWindow<3> (no real cylinder, project to best-fit plane)
//     → LagrangeWindow<3> (conic also invalid; used by blend_curve level)
//
//   The cylinder solver uses 2D Newton on G(a,c)=0, H(a,c)=0 from a 40×40
//   sinh-spaced starting grid (3 axis permutations).  Avoids the Sylvester
//   resultant, which has catastrophic cancellation in double precision.
//
// Usage:
//   auto r = fc::blend_curve<3>(ctrl, times, fc::cylinder_tag{});
//
// Compile:
//   g++ -std=c++17 -O2 -I. -o myapp myapp.cpp

#pragma once
#include "conicblend.hpp"

namespace fc {

struct cylinder_tag {};

namespace detail {

// ── Univariate polynomial utilities ───────────────────────────────────────
static inline double cyl_peval(const double* p, int deg, double t) {
    double v = p[deg];
    for (int i = deg-1; i >= 0; --i) v = v*t + p[i];
    return v;
}

// ── Coefficient structures for F_k(a,b,c,d) = A_k(a,c) + B_k(a,c)·b + D_k(a,c)·d
// Parametrize axis as {y=ax+b, z=cx+d}; differences w.r.t. p[0].
struct CylPDiffs { double X[4],Y[4],Z[4],S[4],P[4],Q[4]; };

static inline void cyl_compB(const CylPDiffs& d, int k, double a, double B[3]) {
    B[0]=2*(d.X[k]*a-d.Y[k]); B[1]=2*d.Z[k]*a; B[2]=-2*d.Y[k];
}
static inline void cyl_compD(const CylPDiffs& d, int k, double a, double D[2]) {
    D[0]=-2*d.Z[k]*(1+a*a); D[1]=2*(d.X[k]+d.Y[k]*a);
}
static inline void cyl_compA(const CylPDiffs& d, int k, double a, double A[3]) {
    double al=d.Y[k]*d.P[k]+d.Z[k]*d.Q[k];
    double be=-d.X[k]*d.P[k]-d.Y[k]*d.S[k];
    double de= d.X[k]*d.S[k]+d.Z[k]*d.Q[k];
    double ga=-d.X[k]*d.Q[k]-d.Z[k]*d.S[k];
    double ze=-d.Y[k]*d.Q[k]-d.Z[k]*d.P[k];
    double ep= d.X[k]*d.S[k]+d.Y[k]*d.P[k];
    A[0]=al+be*a+de*a*a; A[1]=ga+ze*a; A[2]=ep;
}

// Evaluate G(a,c) or H(a,c) = 3×3 minor of [A_r,B_r,D_r] rows
static inline double cyl_eval_minor(const CylPDiffs& d,
                                     double av, double cv,
                                     int r0, int r1, int r2)
{
    double A0[3],B0[3],D0[2], A1[3],B1[3],D1[2], A2[3],B2[3],D2[2];
    cyl_compA(d,r0,av,A0); cyl_compB(d,r0,av,B0); cyl_compD(d,r0,av,D0);
    cyl_compA(d,r1,av,A1); cyl_compB(d,r1,av,B1); cyl_compD(d,r1,av,D1);
    cyl_compA(d,r2,av,A2); cyl_compB(d,r2,av,B2); cyl_compD(d,r2,av,D2);
    double a0=cyl_peval(A0,2,cv), b0=cyl_peval(B0,2,cv), d0=cyl_peval(D0,1,cv);
    double a1=cyl_peval(A1,2,cv), b1=cyl_peval(B1,2,cv), d1=cyl_peval(D1,1,cv);
    double a2=cyl_peval(A2,2,cv), b2=cyl_peval(B2,2,cv), d2=cyl_peval(D2,1,cv);
    return a0*(b1*d2-b2*d1) - b0*(a1*d2-a2*d1) + d0*(a1*b2-a2*b1);
}

// Solve 2×2 linear system for (b,d) given a,c
static inline bool cyl_solve_bd(const CylPDiffs& d, int i, int j,
                                  double a, double c,
                                  double& b_out, double& d_out)
{
    double Ai[3],Bi[3],Di[2],Aj[3],Bj[3],Dj[2];
    cyl_compA(d,i,a,Ai); cyl_compB(d,i,a,Bi); cyl_compD(d,i,a,Di);
    cyl_compA(d,j,a,Aj); cyl_compB(d,j,a,Bj); cyl_compD(d,j,a,Dj);
    double Biv=cyl_peval(Bi,2,c), Div=cyl_peval(Di,1,c);
    double Bjv=cyl_peval(Bj,2,c), Djv=cyl_peval(Dj,1,c);
    double Aiv=cyl_peval(Ai,2,c), Ajv=cyl_peval(Aj,2,c);
    double det=Biv*Djv-Bjv*Div;
    if (std::abs(det)<1e-10) return false;
    b_out=(-Aiv*Djv+Ajv*Div)/det;
    d_out=( Aiv*Bjv-Ajv*Biv)/det;
    return true;
}

struct CylRawSol { double a, c, b, d, rsqr, eq_err; };

// ── 2D Newton solver: for one x-parametrization, find all cylinders ───────
static inline void cyl_solve_xpara(const double pts[5][3],
                                    std::vector<CylRawSol>& out)
{
    CylPDiffs d;
    for (int k=1;k<=4;++k) {
        int j=k-1;
        d.X[j]=pts[k][0]-pts[0][0]; d.Y[j]=pts[k][1]-pts[0][1]; d.Z[j]=pts[k][2]-pts[0][2];
        d.S[j]=pts[k][0]+pts[0][0]; d.P[j]=pts[k][1]+pts[0][1]; d.Q[j]=pts[k][2]+pts[0][2];
    }
    auto G = [&](double av, double cv){ return cyl_eval_minor(d,av,cv,0,1,2); };
    auto H = [&](double av, double cv){ return cyl_eval_minor(d,av,cv,1,2,3); };

    // Back-substitute (a,c) → (b,d,r²), verify, push if accepted
    auto try_ac = [&](double ar, double cr) {
        int bpi=0, bpj=1; double best_det=0;
        const int pairs[6][2]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
        for (auto& pr : pairs) {
            double Bi2[3],Di2[2],Bj2[3],Dj2[2];
            cyl_compB(d,pr[0],ar,Bi2); cyl_compD(d,pr[0],ar,Di2);
            cyl_compB(d,pr[1],ar,Bj2); cyl_compD(d,pr[1],ar,Dj2);
            double dt=std::abs(cyl_peval(Bi2,2,cr)*cyl_peval(Dj2,1,cr)
                              -cyl_peval(Bj2,2,cr)*cyl_peval(Di2,1,cr));
            if (dt>best_det){best_det=dt; bpi=pr[0]; bpj=pr[1];}
        }
        double b_r, d_r;
        if (!cyl_solve_bd(d,bpi,bpj,ar,cr,b_r,d_r)) return;

        double W=1+ar*ar+cr*cr;
        double off[3]={0,b_r,d_r};
        double q0[3]={pts[0][0]-off[0], pts[0][1]-off[1], pts[0][2]-off[2]};
        double qu0=q0[0]+ar*q0[1]+cr*q0[2];
        double rsqr=(q0[0]*q0[0]+q0[1]*q0[1]+q0[2]*q0[2])-qu0*qu0/W;
        if (rsqr<=1e-12) return;

        double max_err=0;
        for (int k=0;k<5;++k) {
            double qk[3]={pts[k][0]-off[0], pts[k][1]-off[1], pts[k][2]-off[2]};
            double quk=qk[0]+ar*qk[1]+cr*qk[2];
            double eq=W*(qk[0]*qk[0]+qk[1]*qk[1]+qk[2]*qk[2]-rsqr)-quk*quk;
            max_err=std::max(max_err,std::abs(eq));
        }
        if (max_err>1e-6*(rsqr*W)+1e-8) return;

        // Dedup: keep best (lowest eq_err) among near-duplicates
        for (auto& s : out) {
            if (std::abs(s.a-ar)<1e-4*(std::abs(s.a)+std::abs(ar)+1) &&
                std::abs(s.c-cr)<1e-4*(std::abs(s.c)+std::abs(cr)+1)) {
                if (max_err<s.eq_err) s={ar,cr,b_r,d_r,rsqr,max_err};
                return;
            }
        }
        out.push_back({ar,cr,b_r,d_r,rsqr,max_err});
    };

    // 2D Newton from each starting point
    auto run_newton = [&](double a0, double c0) {
        double a=a0, c=c0, prev_res=1e300;
        int stall=0;
        for (int iter=0;iter<80;++iter) {
            double Gv=G(a,c), Hv=H(a,c);
            double res=std::abs(Gv)+std::abs(Hv);
            if (res>prev_res*10.){if(++stall>3)break;}else stall=0;
            prev_res=res;
            double es=(res<1e-8)?1e-8:1e-5;
            double ea=es*(1+std::abs(a)), ec=es*(1+std::abs(c));
            double dGda=(G(a+ea,c)-Gv)/ea, dGdc=(G(a,c+ec)-Gv)/ec;
            double dHda=(H(a+ea,c)-Hv)/ea, dHdc=(H(a,c+ec)-Hv)/ec;
            double detJ=dGda*dHdc-dGdc*dHda;
            if (std::abs(detJ)<1e-40) break;
            double da=-(Gv*dHdc-Hv*dGdc)/detJ;
            double dc=-(dGda*Hv-dHda*Gv)/detJ;
            double step=std::sqrt(da*da+dc*dc);
            double maxstep=5.0*(1+std::abs(a)+std::abs(c));
            if (step>maxstep){da*=maxstep/step; dc*=maxstep/step;}
            a+=da; c+=dc;
        }
        try_ac(a,c);
    };

    const int NA=40, NC=40;
    const double ASCALE=200.0, CSCALE=200.0;
    for (int ia=0;ia<NA;++ia) {
        double at=(double)ia/(NA-1)-0.5;
        double a0=ASCALE*std::sinh(5.0*at);
        for (int ic=0;ic<NC;++ic) {
            double ct=(double)ic/(NC-1)-0.5;
            run_newton(a0, CSCALE*std::sinh(5.0*ct));
        }
    }
    for (double av:{0.0,-0.1,0.1,-1.,1.,-5.,5.})
        for (double cv:{0.0,-0.1,0.1,-1.,1.,-5.,5.})
            run_newton(av,cv);
}

// ── Cylinder solution in 3D ────────────────────────────────────────────────
struct CylSol {
    VecN<3> u_hat, v_hat, w_hat, offset;
    double r, phi_span, knot_err;
    double lam_ratio;   // lam_min/lam_max of unrolled 2D scatter; near 0 → geodesic (line)
    bool phi_monotone;  // true iff phi[0..4] is monotone after unwrapping
};

static inline VecN<3> cyl_cross3(VecN<3> const& a, VecN<3> const& b) {
    return VecN<3>{a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}

static inline CylSol cyl_raw_to_sol(const CylRawSol& s, int perm,
                                     const VecN<3> pts5[5])
{
    // Axis direction and offset in permuted coords
    VecN<3> u_p{1, s.a, s.c};
    VecN<3> o_p{0, s.b, s.d};
    u_p = u_p.normalized();

    // Un-permute back to original coords
    // perm=0: (x,y,z) identity
    // perm=1: (y,z,x) → original = (perm[2],perm[0],perm[1])
    // perm=2: (z,x,y) → original = (perm[1],perm[2],perm[0])
    VecN<3> uh, off;
    if (perm==0) { uh=u_p; off=o_p; }
    else if (perm==1) { uh=VecN<3>{u_p[2],u_p[0],u_p[1]}; off=VecN<3>{o_p[2],o_p[0],o_p[1]}; }
    else               { uh=VecN<3>{u_p[1],u_p[2],u_p[0]}; off=VecN<3>{o_p[1],o_p[2],o_p[0]}; }

    // Build in-plane basis (v_hat, w_hat)
    VecN<3> tmp = (std::abs(uh[0])<0.9) ? VecN<3>{1,0,0} : VecN<3>{0,1,0};
    VecN<3> vh = tmp - tmp.dot(uh)*uh; vh = vh.normalized();
    VecN<3> wh = cyl_cross3(uh, vh);

    // Unroll phi and compute phi_span
    double phi[5];
    for (int k=0;k<5;++k) {
        VecN<3> p = pts5[k] - off;
        VecN<3> pp = p - uh*p.dot(uh);
        phi[k] = std::atan2(pp.dot(wh), pp.dot(vh));
    }
    for (int k=1;k<5;++k) {
        double dv=phi[k]-phi[k-1];
        while(dv> M_PI){dv-=2*M_PI; phi[k]-=2*M_PI;}
        while(dv<-M_PI){dv+=2*M_PI; phi[k]+=2*M_PI;}
    }
    double span=std::abs(phi[4]-phi[0]);

    // Phi monotonicity: all consecutive diffs must have the same sign
    bool phi_monotone = true;
    {
        int sign0 = 0;
        for (int k = 1; k < 5; ++k) {
            double dv = phi[k] - phi[k-1];
            if (std::abs(dv) < 1e-9) continue;
            int sgn = (dv > 0) ? 1 : -1;
            if (sign0 == 0) sign0 = sgn;
            else if (sgn != sign0) { phi_monotone = false; break; }
        }
    }

    // Knot error
    double r=std::sqrt(s.rsqr);
    double kerr=0;
    for (int k=0;k<5;++k) {
        VecN<3> q=pts5[k]-off;
        VecN<3> pp=q-uh*q.dot(uh);
        double d=pp.norm()-r; kerr+=d*d;
    }
    kerr=std::sqrt(kerr/5.0);

    // 2D straightness of unrolled points: lam_ratio = lam_min/lam_max.
    // Near 0 → points are near-collinear (curve is a geodesic on this cylinder,
    // e.g. helix on its own cylinder).  Used to prefer the natural cylinder.
    double lam_ratio = 1.0;
    {
        double mx = 0, mz = 0;
        for (int k = 0; k < 5; ++k) {
            VecN<3> q = pts5[k] - off;
            mx += r * phi[k];
            mz += q.dot(uh);
        }
        mx /= 5; mz /= 5;
        double sxx=0,sxz=0,szz=0;
        for (int k=0;k<5;++k) {
            VecN<3> q = pts5[k] - off;
            double dx = r*phi[k] - mx, dz = q.dot(uh) - mz;
            sxx+=dx*dx; sxz+=dx*dz; szz+=dz*dz;
        }
        double tr = sxx + szz;
        double dsc = (sxx-szz)*(sxx-szz) + 4*sxz*sxz;
        double lm = 0.5*(tr - std::sqrt(dsc));
        double lM = 0.5*(tr + std::sqrt(dsc));
        lam_ratio = (lM > 1e-30) ? lm / lM : 0.0;
    }

    return CylSol{uh, vh, wh, off, r, span, kerr, lam_ratio, phi_monotone};
}

static inline bool cyl_same(const CylSol& a, const CylSol& b) {
    if (std::abs(a.r-b.r)>0.02*(a.r+b.r+1e-10)) return false;
    return std::abs(a.u_hat.dot(b.u_hat))>0.9998;
}

// Find all real right circular cylinders through 5 points.
// Returns solutions sorted: monotone-phi first, then by ascending phi_span.
// Non-monotone cylinders are still returned (at the end) as a last resort.
static inline std::vector<CylSol> cyl_solve(const VecN<3> pts5[5])
{
    static const int perm_map[3][3]={{0,1,2},{1,2,0},{2,0,1}};
    std::vector<CylSol> all;
    for (int perm=0;perm<3;++perm) {
        double p[5][3];
        for (int i=0;i<5;++i)
            for (int dd=0;dd<3;++dd)
                p[i][dd]=pts5[i][perm_map[perm][dd]];
        std::vector<CylRawSol> raw;
        cyl_solve_xpara(p, raw);
        for (auto& rs : raw) {
            CylSol cs=cyl_raw_to_sol(rs,perm,pts5);
            bool dup=false;
            for (auto& ex : all) if(cyl_same(ex,cs)){dup=true;break;}
            if (!dup) all.push_back(cs);
        }
    }
    // Sort: monotone-phi first.  A non-monotone phi sequence means the control
    // points' angular projections fold back as t increases — the cylinder is a
    // non-injective frame for this arc.  Within monotone, geodesic cylinders
    // (lam_ratio < 1e-4 → near-collinear unrolled points, e.g. helix on its
    // own cylinder) come first because they give exact Lagrange interpolation.
    // Among non-geodesic cylinders, sort by phi_span (smaller = more compact).
    std::sort(all.begin(),all.end(),[](const CylSol&a,const CylSol&b){
        if (a.phi_monotone != b.phi_monotone) return a.phi_monotone > b.phi_monotone;
        bool ag = (a.lam_ratio < 1e-4), bg = (b.lam_ratio < 1e-4);
        if (ag != bg) return ag > bg;
        if (ag && bg) return a.lam_ratio < b.lam_ratio;
        return a.phi_span < b.phi_span;
    });
    return all;
}

} // namespace detail

// ── CylinderWindow<3> ─────────────────────────────────────────────────────
//
// Fits a right circular cylinder through 5 control points (smallest Δφ),
// unrolls to (r·φ, z), and applies a 2D ConicWindow there.
// Falls back to ConicWindow<3> if no real cylinder is found.
// valid() is false only if both paths fail → blend_curve uses LagrangeWindow.

template <int Dim>
class CylinderWindow;   // defined only for Dim=3

template <>
class CylinderWindow<3> {
    bool valid_       = false;
    bool use_cylinder_= false;

    // Cylinder path: geometry for re-rolling + inner 2D window
    // When the unrolled points are near-collinear (e.g. helix), ConicWindow<2>
    // fails its collinearity check; we fall back to LagrangeWindow<2> in the
    // cylinder frame, which reproduces geodesics (e.g. helix) exactly.
    VecN<3> u_hat_{}, v_hat_{}, w_hat_{}, offset_{};
    double  r_ = 0.0;
    using Inner2D = std::variant<ConicWindow<2>, LagrangeWindow<2>>;
    Inner2D inner_2d_;

    // Fallback path: plane-projected conic (when no real cylinder found)
    ConicWindow<3> conic3d_;

public:
    CylinderWindow() = default;

    CylinderWindow(VecN<3> const& p0, VecN<3> const& p1, VecN<3> const& p2,
                   VecN<3> const& p3, VecN<3> const& p4,
                   double t0, double t1, double t2, double t3, double t4)
    {
        VecN<3> pts5[5]={p0,p1,p2,p3,p4};

        // ── Try cylinder path ──────────────────────────────────────────────
        auto cyls = detail::cyl_solve(pts5);
        if (!cyls.empty()) {
            auto& cs = cyls[0];  // smallest phi_span

            u_hat_ = cs.u_hat; v_hat_ = cs.v_hat;
            w_hat_ = cs.w_hat; offset_ = cs.offset; r_ = cs.r;

            // Unroll 5 control points to (r·φ, z)
            VecN<2> upts[5];
            double phi[5];
            for (int k=0;k<5;++k) {
                VecN<3> p = pts5[k] - offset_;
                VecN<3> pp = p - u_hat_*p.dot(u_hat_);
                phi[k] = std::atan2(pp.dot(w_hat_), pp.dot(v_hat_));
            }
            for (int k=1;k<5;++k) {    // unwrap
                double dv=phi[k]-phi[k-1];
                while(dv> M_PI){dv-=2*M_PI; phi[k]-=2*M_PI;}
                while(dv<-M_PI){dv+=2*M_PI; phi[k]+=2*M_PI;}
            }
            for (int k=0;k<5;++k) {
                double z = (pts5[k]-offset_).dot(u_hat_);
                upts[k] = VecN<2>{r_*phi[k], z};
            }

            // Try 2D conic; if unrolled points are near-collinear (e.g. helix
            // is a geodesic → straight line after unrolling), fall back to
            // Lagrange<2> which reproduces a line exactly via polynomial interp.
            ConicWindow<2> cw2(upts[0],upts[1],upts[2],upts[3],upts[4],
                               t0,t1,t2,t3,t4);
            if (cw2.valid()) {
                inner_2d_ = std::move(cw2);
            } else {
                inner_2d_ = LagrangeWindow<2>(upts[0],upts[1],upts[2],upts[3],upts[4],
                                              t0,t1,t2,t3,t4);
            }
            use_cylinder_ = true;
            valid_ = true;
            return;
        }

        // ── Fallback: ConicWindow<3> on best-fit plane ────────────────────
        conic3d_ = ConicWindow<3>(p0,p1,p2,p3,p4, t0,t1,t2,t3,t4);
        if (conic3d_.valid()) {
            valid_ = true;
        }
    }

    bool valid()         const { return valid_; }
    bool used_cylinder() const { return use_cylinder_; }

    VecN<3> operator()(double t) const {
        if (use_cylinder_) {
            VecN<2> uv = std::visit([t](auto const& w) { return w(t); }, inner_2d_);
            double phi = uv[0] / r_;
            double z   = uv[1];
            return offset_
                 + u_hat_ * z
                 + v_hat_ * (r_ * std::cos(phi))
                 + w_hat_ * (r_ * std::sin(phi));
        }
        return conic3d_(t);
    }
};

// ── blend_curve(..., cylinder_tag{}) ──────────────────────────────────────
//
// Builds CylinderWindow<3> per 5-point window.  Falls back to LagrangeWindow
// for windows where both cylinder and planar conic fail.
//
// used_cylinder (optional): set to true if at least one cylinder window was used.

template <int Dim>
inline BlendResultND<Dim> blend_curve(
#if __cplusplus >= 202002L
    std::span<VecN<Dim> const> ctrl,
    std::span<double const>    times,
#else
    std::vector<VecN<Dim>> const& ctrl,
    std::vector<double>    const& times,
#endif
    cylinder_tag,
    int   pts_per_seg   = 60,
    int   smooth_N      = 2,
    bool* used_cylinder = nullptr)
{
    static_assert(Dim == 3, "cylinder_tag blend_curve is only defined for Dim=3");

    int n = static_cast<int>(ctrl.size());
    if (n != static_cast<int>(times.size()))
        FC_THROW("fc::blend_curve(cylinder): ctrl and times must have equal length");
    if (n < 6)
        FC_THROW("fc::blend_curve(cylinder): need at least 6 control points");
    if (pts_per_seg < 2)
        FC_THROW("fc::blend_curve(cylinder): pts_per_seg must be >= 2");

    using Win = std::variant<CylinderWindow<3>, LagrangeWindow<3>>;
    std::vector<Win> wins;
    wins.reserve(n-4);
    int n_cyl = 0;

    for (int i=0;i<n-4;++i) {
        CylinderWindow<3> cw(ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4],
                             times[i],times[i+1],times[i+2],times[i+3],times[i+4]);
        if (cw.valid()) {
            if (cw.used_cylinder()) ++n_cyl;
            wins.emplace_back(std::move(cw));
        } else {
            wins.emplace_back(LagrangeWindow<3>(
                ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4],
                times[i],times[i+1],times[i+2],times[i+3],times[i+4]));
        }
    }
    if (used_cylinder) *used_cylinder = (n_cyl > 0);

    auto eval_win = [](Win const& w, double t) -> VecN<3> {
        return std::visit([t](auto const& ww) -> VecN<3> { return ww(t); }, w);
    };

    int n_segs = n-5;
    BlendResultND<3> out;
    out.pts.reserve(n_segs*pts_per_seg+1);
    out.times.reserve(n_segs*pts_per_seg+1);

    for (int j=2;j<=n-4;++j) {
        Win const& wA=wins[j-2], &wB=wins[j-1];
        double t_lo=times[j], t_hi=times[j+1];
        int k_max=(j<n-4)?pts_per_seg-1:pts_per_seg;
        for (int k=0;k<=k_max;++k) {
            double s=(double)k/pts_per_seg;
            double t=t_lo+s*(t_hi-t_lo);
            double w=smoothstep(s,smooth_N);
            VecN<3> A=eval_win(wA,t), B=eval_win(wB,t);
            out.pts.push_back(A*(1.0-w)+B*w);
            out.times.push_back(t);
        }
    }
    return out;
}

// 3D non-template convenience overload
inline BlendResult blend_curve(
#if __cplusplus >= 202002L
    std::span<Vec3 const>   ctrl,
    std::span<double const> times,
#else
    std::vector<Vec3>   const& ctrl,
    std::vector<double> const& times,
#endif
    cylinder_tag,
    int   pts_per_seg   = 60,
    int   smooth_N      = 2,
    bool* used_cylinder = nullptr)
{
    return blend_curve<3>(ctrl, times, cylinder_tag{},
                          pts_per_seg, smooth_N, used_cylinder);
}

} // namespace fc
