// conicblend_cylinder.hpp — C^N 3D curve interpolation via cylinder unrolling
// Header-only C++17 library.  Requires conicblend.hpp.
//
// Architecture:
//   Each 5-point window tries all real cylinders through the 5 control points
//   in sorted order (geodesic first, then phi_span ascending).  For each,
//   ConicWindow<2> is attempted on the unrolled (r·φ, z) control points:
//     - geodesic (lam_ratio < 1e-4): unrolled pts are near-collinear → ConicWindow<2>
//       uses its internal line mode (degree-4 Lagrange) — machine-precision for helices.
//     - non-geodesic: ConicWindow<2> uses the standard conic rational map.
//   If ConicWindow<2> rejects (degenerate conic, non-monotone arc-angle), the next
//   cylinder is tried.  The re-rolled result is the final 3D curve for this window.
//
//   Fallback chain per window:
//     Exact-fit cylinder + ConicWindow<2>     (5 pts exactly on cylinder)
//     → Best-fit cylinder + ConicWindow<2>    (only if max surface residual < 5% of
//                                              window scale; grid search on S²)
//     → LagrangeWindow<3>                     (both fail; used by blend_curve level)
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
    // Sort: geodesic cylinders (lam_ratio < 1e-4 → unrolled points near-collinear,
    // e.g. helix on its own cylinder) first — ConicWindow<2> uses line mode for
    // these, giving machine-precision interpolation.
    // Among non-geodesic cylinders, sort by phi_span ascending (smaller = more
    // compact = more numerically stable conic fit).
    //
    // phi_monotone is NOT used as a sort key here.  It correctly guards
    // LagrangeWindow<2> (a non-monotone r·φ(t) signal makes Lagrange oscillate),
    // but it is the wrong gate for ConicWindow<2>, which only requires that the
    // conic arc-angle from P₀ is monotone — a strictly weaker condition that
    // ConicWindow<2> checks internally.  CylinderWindow iterates through this
    // sorted list and lets ConicWindow<2>'s own validity test decide.
    std::sort(all.begin(),all.end(),[](const CylSol&a,const CylSol&b){
        bool ag = (a.lam_ratio < 1e-4), bg = (b.lam_ratio < 1e-4);
        if (ag != bg) return ag > bg;
        if (ag && bg) return a.lam_ratio < b.lam_ratio;
        return a.phi_span < b.phi_span;
    });
    return all;
}

// ── Best-fit right circular cylinder ─────────────────────────────────────
// Minimise Σ(dist(pᵢ,axis) − r)² with r = mean dist (eliminated analytically).
// Axis direction u searched over S² on a 20×20 grid; for each u the optimal
// centre and r are found via Coope's algebraic circle fit (LS on the linear
// system a·xᵢ + b·yᵢ + c = xᵢ²+yᵢ²) applied to the 2D projections ⊥ u.
//
// Returns the best (u, centre, r) found, plus max_surf_res = normalised max
// surface residual max(|dist(pᵢ,axis)−r|) / win_scale.  Caller gates on this.

struct CylBFSol {
    VecN<3> u_hat, v_hat, w_hat, offset;
    double  r, max_surf_res;  // max_surf_res dimensionless; 1e30 if degenerate
};

static inline CylBFSol cyl_best_fit(const VecN<3> pts5[5], double win_scale)
{
    CylBFSol best;
    best.max_surf_res = 1e30;

    auto try_u = [&](VecN<3> u) {
        double un = u.norm();
        if (un < 1e-10) return;
        u = u * (1.0 / un);

        // Perpendicular basis
        VecN<3> e1 = (std::abs(u[0]) < 0.9) ? VecN<3>{1,0,0} : VecN<3>{0,1,0};
        e1 = e1 - u * e1.dot(u);
        double e1n = e1.norm();
        if (e1n < 1e-10) return;
        e1 = e1 * (1.0 / e1n);
        VecN<3> e2 = cyl_cross3(u, e1);

        // Project 5 points to 2D plane ⊥ u
        double px[5], py[5];
        for (int i = 0; i < 5; ++i) {
            px[i] = pts5[i].dot(e1);
            py[i] = pts5[i].dot(e2);
        }

        // Coope's algebraic circle fit: solve a·xᵢ + b·yᵢ + c = zᵢ (LS) where
        // zᵢ = xᵢ²+yᵢ².  Center = (a/2, b/2), r² = c + cx²+cy².
        // This minimises Σ(algebraic distance)² and gives the correct circle for
        // exact data; it does NOT minimise Σ(dᵢ−r̄)² but is an excellent proxy
        // and far cheaper than a geometric iteration.
        double M[3][3]={}, rhs[3]={};
        for (int i = 0; i < 5; ++i) {
            double x=px[i], y=py[i], z=x*x+y*y;
            M[0][0]+=x*x; M[0][1]+=x*y; M[0][2]+=x;
            M[1][0]+=x*y; M[1][1]+=y*y; M[1][2]+=y;
            M[2][0]+=x;   M[2][1]+=y;   M[2][2]+=1;
            rhs[0]+=z*x;  rhs[1]+=z*y;  rhs[2]+=z;
        }
        // 3×3 Gaussian elimination (no pivoting; well-conditioned for spread pts)
        double sol[3]={rhs[0],rhs[1],rhs[2]};
        for (int col = 0; col < 3; ++col) {
            if (std::abs(M[col][col]) < 1e-30) return;   // degenerate
            for (int row = col+1; row < 3; ++row) {
                double f = M[row][col] / M[col][col];
                for (int k = col; k < 3; ++k) M[row][k] -= f * M[col][k];
                sol[row] -= f * sol[col];
            }
        }
        // Back substitution
        for (int row = 2; row >= 0; --row) {
            for (int k = row+1; k < 3; ++k) sol[row] -= M[row][k] * sol[k];
            sol[row] /= M[row][row];
        }
        double cx = sol[0]*0.5, cy = sol[1]*0.5;
        double rsq = sol[2] + cx*cx + cy*cy;
        if (rsq <= 0.0) return;
        double r = std::sqrt(rsq);
        if (r < 1e-10) return;

        // Max surface residual (normalised)
        double maxres = 0;
        for (int i = 0; i < 5; ++i) {
            double d = std::sqrt((px[i]-cx)*(px[i]-cx) + (py[i]-cy)*(py[i]-cy));
            double res = std::abs(d - r);
            if (res > maxres) maxres = res;
        }
        double rel = (win_scale > 0) ? maxres / win_scale : maxres;
        if (rel < best.max_surf_res) {
            best.max_surf_res = rel;
            best.r      = r;
            best.u_hat  = u;
            best.v_hat  = e1;
            best.w_hat  = e2;
            best.offset = e1 * cx + e2 * cy;   // point on axis with zero u-component
        }
    };

    // Grid search on S²
    const int NT = 20, NP = 20;
    for (int it = 0; it < NT; ++it) {
        double theta = M_PI * (it + 0.5) / NT;
        double st = std::sin(theta), ct = std::cos(theta);
        for (int ip = 0; ip < NP; ++ip) {
            double phi = 2.0 * M_PI * ip / NP;
            try_u(VecN<3>{st*std::cos(phi), st*std::sin(phi), ct});
        }
    }
    // Axis-aligned extras
    for (auto ax : std::initializer_list<VecN<3>>{{1,0,0},{0,1,0},{0,0,1}})
        try_u(ax);

    return best;
}

} // namespace detail

// ── CylinderWindow<3> ─────────────────────────────────────────────────────
//
// Fits a cylinder through 5 control points, unrolls to (r·φ, z), and applies
// ConicWindow<2> there.  Fallback chain:
//   1. Exact-fit cylinder (Newton solver) + ConicWindow<2>
//   2. Best-fit cylinder (grid search on S²) + ConicWindow<2>,
//      only if max surface residual < 5% of window scale
//   3. Caller falls back to LagrangeWindow<3> (valid() == false)

template <int Dim>
class CylinderWindow;   // defined only for Dim=3

template <>
class CylinderWindow<3> {
    bool valid_        = false;
    bool use_best_fit_ = false;   // true → level 2 (best-fit) was used

    VecN<3> u_hat_{}, v_hat_{}, w_hat_{}, offset_{};
    double  r_ = 0.0;
    using Inner2D = ConicWindow<2>;
    Inner2D inner_2d_;

    // Shared helper: unroll pts5 onto a cylinder (u_hat, v_hat, w_hat, offset, r)
    // and attempt ConicWindow<2>.  Fills members and returns true on success.
    bool try_cylinder_(const VecN<3> pts5[5],
                       VecN<3> const& u, VecN<3> const& v, VecN<3> const& w,
                       VecN<3> const& off, double r,
                       double t0, double t1, double t2, double t3, double t4)
    {
        VecN<2> upts[5];
        double phi[5];
        for (int k = 0; k < 5; ++k) {
            VecN<3> p  = pts5[k] - off;
            VecN<3> pp = p - u * p.dot(u);
            phi[k] = std::atan2(pp.dot(w), pp.dot(v));
        }
        for (int k = 1; k < 5; ++k) {
            double dv = phi[k] - phi[k-1];
            while (dv >  M_PI) { dv -= 2*M_PI; phi[k] -= 2*M_PI; }
            while (dv < -M_PI) { dv += 2*M_PI; phi[k] += 2*M_PI; }
        }
        for (int k = 0; k < 5; ++k) {
            double z = (pts5[k] - off).dot(u);
            upts[k] = VecN<2>{r * phi[k], z};
        }
        Inner2D win(upts[0],upts[1],upts[2],upts[3],upts[4],
                    t0,t1,t2,t3,t4);
        if (!win.valid()) return false;
        u_hat_ = u; v_hat_ = v; w_hat_ = w; offset_ = off; r_ = r;
        inner_2d_ = std::move(win);
        return true;
    }

public:
    CylinderWindow() = default;

    CylinderWindow(VecN<3> const& p0, VecN<3> const& p1, VecN<3> const& p2,
                   VecN<3> const& p3, VecN<3> const& p4,
                   double t0, double t1, double t2, double t3, double t4)
    {
        VecN<3> pts5[5] = {p0,p1,p2,p3,p4};

        // ── Level 1: exact-fit cylinders (Newton solver) ─────────────────
        // Sorted: geodesic first, then phi_span ascending.
        // ConicWindow<2> handles both line-mode (geodesic) and conic cases.
        auto cyls = detail::cyl_solve(pts5);
        for (auto& cs : cyls) {
            if (try_cylinder_(pts5, cs.u_hat, cs.v_hat, cs.w_hat, cs.offset, cs.r,
                              t0,t1,t2,t3,t4)) {
                valid_ = true;
                return;
            }
        }

        // ── Level 2: best-fit cylinder ────────────────────────────────────
        // Gate: only proceed if max surface residual < 5% of window scale.
        // Invariant: residual/win_scale is dimensionless, rigid-motion invariant.
        double win_scale = 0.0;
        for (int i = 0; i < 5; ++i)
            for (int j = i+1; j < 5; ++j)
                win_scale = std::max(win_scale, (pts5[i]-pts5[j]).norm());

        auto bf = detail::cyl_best_fit(pts5, win_scale);
        if (bf.max_surf_res < 5e-2) {
            if (try_cylinder_(pts5, bf.u_hat, bf.v_hat, bf.w_hat, bf.offset, bf.r,
                              t0,t1,t2,t3,t4)) {
                use_best_fit_ = true;
                valid_        = true;
                return;
            }
        }
        // valid_ stays false → blend_curve falls back to LagrangeWindow<3>
    }

    bool valid()          const { return valid_; }
    bool used_cylinder()  const { return valid_; }   // always true when valid
    bool used_best_fit()  const { return use_best_fit_; }

    VecN<3> operator()(double t) const {
        VecN<2> uv  = inner_2d_(t);
        double  phi = uv[0] / r_;
        double  z   = uv[1];
        return offset_
             + u_hat_ * z
             + v_hat_ * (r_ * std::cos(phi))
             + w_hat_ * (r_ * std::sin(phi));
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
