// diag_parabola2.cpp — step-by-step trace of ConicWindow for parabola windows 3-6
//
// Compile: g++ -std=c++17 -O2 -I. -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 -o diag_parabola2 diag_parabola2.cpp

#include "conicblend.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

using namespace fc;
using namespace fc::detail;
using V2 = VecN<2>;

static constexpr double MY_PI = 3.14159265358979323846;

static std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i) v[i] = a + (b - a) * i / (n - 1);
    return v;
}

static void trace_window(int win_idx, double t0, double t1, double t2, double t3, double t4,
                          V2 p0, V2 p1, V2 p2, V2 p3, V2 p4)
{
    printf("\n=== Window %d: t=[%.4f, %.4f, %.4f, %.4f, %.4f] ===\n",
           win_idx, t0, t1, t2, t3, t4);

    V2 pts5[5] = {p0,p1,p2,p3,p4};
    double ts[5] = {t0,t1,t2,t3,t4};

    // Step 1: project to best-fit plane (trivial for 2D)
    PlaneND<2> plane = best_fit_plane<2>(pts5);
    double (*pts2d)[2] = plane.pts2d;
    printf("  pts2d: ");
    for (int i=0;i<5;i++) printf("(%.4f,%.4f) ", pts2d[i][0], pts2d[i][1]);
    printf("\n");

    // Step 1b: collinearity guard
    {
        double mx=0,my=0;
        for(int i=0;i<5;++i){mx+=pts2d[i][0];my+=pts2d[i][1];}
        mx/=5; my/=5;
        double sxx=0,sxy=0,syy=0;
        for(int i=0;i<5;++i){
            double dx=pts2d[i][0]-mx, dy=pts2d[i][1]-my;
            sxx+=dx*dx; sxy+=dx*dy; syy+=dy*dy;
        }
        double tr=sxx+syy, dsc=(sxx-syy)*(sxx-syy)+4*sxy*sxy;
        double lam_min = 0.5*(tr - std::sqrt(dsc));
        double lam_max = 0.5*(tr + std::sqrt(dsc));
        printf("  collinear: lam_min=%.3e lam_max=%.3e ratio=%.3e %s\n",
               lam_min, lam_max, lam_min/lam_max,
               (lam_max < 1e-30 || lam_min < 1e-4 * lam_max) ? "FAIL" : "ok");
    }

    // Step 2: conic fit
    double coeffs[6];
    fit_conic_2d(pts2d, coeffs);
    double A=coeffs[0], B=coeffs[1], C=coeffs[2], D=coeffs[3], E=coeffs[4], F=coeffs[5];
    printf("  conic [A,B,C,D,E,F] = [%.4e, %.4e, %.4e, %.4e, %.4e, %.4e]\n",
           A, B, C, D, E, F);

    // Step 3: principal frame
    Eig2x2 epf = eigh2(A, B*0.5, C);
    int i0 = (std::abs(epf.val[0]) <= std::abs(epf.val[1])) ? 0 : 1;
    int i1 = 1 - i0;
    double lam0 = epf.val[i0], lam1 = epf.val[i1];
    printf("  eigenvalues: lam0=%.6e lam1=%.6e  product=%.3e\n", lam0, lam1, lam0*lam1);

    double vx0=epf.vec[i0][0], vy0=epf.vec[i0][1];
    double vx1=epf.vec[i1][0], vy1=epf.vec[i1][1];
    // sign convention
    double phys_x = vx0*plane.e1[0]+vy0*plane.e2[0];
    double phys_y = vx0*plane.e1[1]+vy0*plane.e2[1];
    if (phys_x < 0.0 || (std::abs(phys_x)<1e-12 && phys_y<0.0)) { vx0=-vx0; vy0=-vy0; }
    if (vx0*vy1-vy0*vx1 < 0.0) { vx1=-vx1; vy1=-vy1; }

    double pts_p[5][2];
    for (int i=0;i<5;++i) {
        pts_p[i][0] = pts2d[i][0]*vx0 + pts2d[i][1]*vy0;
        pts_p[i][1] = pts2d[i][0]*vx1 + pts2d[i][1]*vy1;
    }
    double A_p=lam0, C_p=lam1;
    double DE[2]={D,E};
    double D_p=DE[0]*vx0+DE[1]*vy0, E_p=DE[0]*vx1+DE[1]*vy1;
    printf("  principal: A_p=%.6e C_p=%.6e D_p=%.6e E_p=%.6e\n", A_p,C_p,D_p,E_p);

    // Step 4: swap
    double A_e=A_p, C_e=C_p, D_e=D_p, E_e=E_p;
    double pts_e[5][2];
    for (int i=0;i<5;++i) { pts_e[i][0]=pts_p[i][0]; pts_e[i][1]=pts_p[i][1]; }
    bool swapped = false;
    constexpr double EPS_SWAP=1e-9;
    if (std::abs(A_e) < EPS_SWAP * std::max(std::abs(C_e),1.0)) {
        std::swap(A_e,C_e); std::swap(D_e,E_e);
        for (int i=0;i<5;++i) std::swap(pts_e[i][0],pts_e[i][1]);
        swapped=true;
    }
    printf("  after swap(%s): A_e=%.6e C_e=%.6e D_e=%.6e E_e=%.6e\n",
           swapped?"yes":"no", A_e,C_e,D_e,E_e);
    printf("  pts_e: ");
    for (int i=0;i<5;i++) printf("(%.4f,%.4f) ", pts_e[i][0], pts_e[i][1]);
    printf("\n");

    // Step 4b: cross-branch
    bool is_cross_branch = false;
    if (lam0*lam1 < 0.0) {
        printf("  -> hyperbola branch detected (lam0*lam1<0)\n");
        double det4 = 4.0*A*C - B*B;
        double cx = (-2.0*C*D+B*E)/det4, cy = (-2.0*A*E+B*D)/det4;
        printf("     center=(%.4e, %.4e)  det4=%.4e\n", cx, cy, det4);
        double tx=(lam0*F<0.0)?vx0:vx1, ty=(lam0*F<0.0)?vy0:vy1;
        for (int i=0; i<4 && !is_cross_branch; ++i) {
            double s0=(pts2d[i][0]-cx)*tx+(pts2d[i][1]-cy)*ty;
            double s1=(pts2d[i+1][0]-cx)*tx+(pts2d[i+1][1]-cy)*ty;
            if (s0*s1<0.0) { is_cross_branch=true; printf("     sign change at i=%d\n",i); }
        }
    }
    printf("  is_cross_branch=%s\n", is_cross_branch?"YES":"no");

    // Step 5a: vertex
    constexpr double EPS_VERTEX=1e-9, EPS_COEFF=1e-12;
    bool has_vertex=false;
    double x0_vert=0, y0_vert=0;
    if (std::abs(A_e) >= EPS_VERTEX * std::max(std::abs(C_e),1.0)) {
        double x_v = -D_e/(2.0*A_e);
        double cv_const = A_e*x_v*x_v + D_e*x_v + F;
        printf("  vertex: x_v=%.6f cv_const=%.6e |C_e|=%.6e |E_e|=%.6e\n",
               x_v, cv_const, std::abs(C_e), std::abs(E_e));
        if (std::abs(C_e) < EPS_VERTEX) {
            if (std::abs(E_e) >= EPS_VERTEX) {
                double y_v = -cv_const / E_e;
                double Mv = std::abs(2.0*C_e*y_v + E_e);
                printf("     parabola branch: y_v=%.6f |M_v|=%.6e\n", y_v, Mv);
                if (Mv >= EPS_COEFF) { x0_vert=x_v; y0_vert=y_v; has_vertex=true; }
            }
        } else {
            double disc = E_e*E_e - 4.0*C_e*cv_const;
            printf("     ellipse branch: disc=%.6e\n", disc);
            if (disc >= 0.0) {
                double sq=std::sqrt(disc);
                double y1=(-E_e+sq)/(2.0*C_e), y2=(-E_e-sq)/(2.0*C_e);
                double cy_mean=0;
                for (int i=0;i<5;++i) cy_mean+=pts_e[i][1]; cy_mean/=5;
                double y_v=(std::abs(y1-cy_mean)<=std::abs(y2-cy_mean))?y1:y2;
                double Mv=std::abs(2.0*C_e*y_v+E_e);
                printf("     y1=%.4f y2=%.4f y_v=%.4f |M_v|=%.4e\n",y1,y2,y_v,Mv);
                if (Mv>=EPS_COEFF) { x0_vert=x_v; y0_vert=y_v; has_vertex=true; }
            }
        }
    } else {
        printf("  vertex: A_e too small (%.6e), skipping\n", std::abs(A_e));
    }
    printf("  has_vertex=%s (%g, %g)\n", has_vertex?"YES":"no", x0_vert, y0_vert);

    // Step 5b: phi_mono lambda (replicated)
    constexpr double EPS_DX=1e-9;
    auto do_phi = [&](double x0, double y0, double Le, double Me, int k_self) -> bool {
        double sv[5], pv[5];
        for (int i=0;i<5;++i) {
            if (i==k_self) sv[i] = -Le/Me;
            else {
                double dx=pts_e[i][0]-x0, dy=pts_e[i][1]-y0;
                if (std::abs(dx) >= EPS_DX*(std::abs(dy)+1.0)) sv[i]=dy/dx;
                else sv[i]=(dy>=0?1.0:-1.0)*1e15;
            }
        }
        for (int i=0;i<5;++i) pv[i]=2.0*std::atan(sv[i]);
        for (int i=0;i<4;++i) {
            double d=pv[i+1]-pv[i];
            while(d>MY_PI){d-=2*MY_PI;pv[i+1]-=2*MY_PI;}
            while(d<-MY_PI){d+=2*MY_PI;pv[i+1]+=2*MY_PI;}
        }
        printf("    sv="); for(int i=0;i<5;i++) printf("%.4f ",sv[i]); printf("\n");
        printf("    pv="); for(int i=0;i<5;i++) printf("%.4f ",pv[i]); printf("\n");
        bool mu=true,md=true;
        for(int i=0;i<4;++i){
            if(pv[i+1]<=pv[i]) mu=false;
            if(pv[i+1]>=pv[i]) md=false;
        }
        return mu||md;
    };

    // Step 5c: try vertex
    bool found=false;
    if (has_vertex) {
        double M_v=2.0*C_e*y0_vert+E_e;
        printf("  try vertex P0=(%.4f,%.4f) M_v=%.4f:\n", x0_vert, y0_vert, M_v);
        bool mono = do_phi(x0_vert, y0_vert, 0.0, M_v, -1);
        printf("    phi_mono=%s\n", mono?"YES":"no");
        if (mono) found=true;
    }

    // Step 5d: fallback argmin
    double L_all[5], M_all[5];
    for (int k=0;k<5;++k) {
        L_all[k]=2.0*A_e*pts_e[k][0]+D_e;
        M_all[k]=2.0*C_e*pts_e[k][1]+E_e;
    }
    if (!found) {
        double best=1e30; int k_best=-1;
        for (int k=1;k<=3;++k) {
            double Mk=M_all[k];
            if (std::abs(Mk)<EPS_COEFF) continue;
            double r=std::abs(L_all[k]/Mk);
            if (r<best){best=r;k_best=k;}
        }
        if (k_best>=0) {
            double Lk=L_all[k_best], Mk=M_all[k_best];
            printf("  try argmin k=%d P0=(%.4f,%.4f) L/M=%.4f:\n",
                   k_best, pts_e[k_best][0], pts_e[k_best][1], Lk/Mk);
            bool mono=do_phi(pts_e[k_best][0], pts_e[k_best][1], Lk, Mk, k_best);
            printf("    phi_mono=%s\n", mono?"YES":"no");
            if (mono) found=true;
        }
    }

    // Cross-branch gate
    if (found && is_cross_branch) {
        printf("  cross-branch gate: REJECTED (allow_cross_branch=false)\n");
        found=false;
    }

    ConicWindow<2> w(p0,p1,p2,p3,p4,t0,t1,t2,t3,t4);
    printf("  => found=%s  ConicWindow::valid=%s\n",
           found?"yes":"NO", w.valid()?"YES":"FAIL");
}

int main() {
    int n = 14;
    auto ts = linspace(-2.0, 2.0, n);
    std::vector<V2> ctrl(n);
    for (int i=0;i<n;++i) ctrl[i]=V2(ts[i],ts[i]*ts[i]);

    for (int i : {3, 4, 5, 6}) {
        trace_window(i, ts[i],ts[i+1],ts[i+2],ts[i+3],ts[i+4],
                     ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4]);
    }
    return 0;
}
