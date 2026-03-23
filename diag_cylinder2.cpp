// diag_cylinder2.cpp — Lichtblau polynomial cylinder solver (C++ port of cylinder_solver.wls)
//
// Finds all real right circular cylinders through 5 3D points.
// Algorithm (Lichtblau ICMS 2006):
//   Axis parametrized as {y=ax+b, z=cx+d}.  For each of 3 permutations of (x,y,z) we
//   apply the "x-parametrization": axis direction (1,a,c), offset (0,b,d).
//   The 5 distance equations eliminate r² via differences → 4 equations F_k(a,b,c,d)=0,
//   each linear in (b,d).  Form the 3×3 minors
//       G = det([F1;F2;F3]) and H = det([F2;F3;F4])  (bivariate polys in a,c)
//   and solve G(a,c)=0, H(a,c)=0 simultaneously via 2D Newton from a dense sinh-spaced
//   starting grid.  Back-substitute for (b,d,r²) and verify vs all 5 equations.
//
//   Note: the classical Sylvester resultant approach (Res_c(G,H) as a univariate poly in a)
//   suffers catastrophic cancellation in double precision; 2D Newton avoids this entirely.
//
// Compile:
//   g++ -std=c++17 -O2 -I. -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 \
//       -o diag_cylinder2 diag_cylinder2.cpp

#include "conicblend.hpp"
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <functional>
#include <cstring>

using V3 = fc::VecN<3>;

static V3 cross3(V3 const& a, V3 const& b) {
    return V3(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}
static double vnorm2(V3 const& v){return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];}
static double vdot(V3 const& a, V3 const& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static double vnorm(V3 const& v){return std::sqrt(vnorm2(v));}

// ── Univariate polynomial utilities ───────────────────────────────────────
// Poly: p[0]+p[1]*t+...  size = degree+1


static double peval(const double* p, int deg, double t) {
    double v = p[deg];
    for (int i=deg-1;i>=0;--i) v=v*t+p[i];
    return v;
}

// ── A_k, B_k, D_k polynomials in c (for fixed a) ─────────────────────────
// Using differences w.r.t. p0: X[k]=p[k+1][0]-p[0][0], Y[k]=p[k+1][1]-p[0][1], etc.
// (k is 0-indexed: k=0..3 corresponds to original differences 1..4)
//
// B_k = [2*(Xk*a-Yk), 2*Zk*a, -2*Yk]        (degree 2)
// D_k = [-2*Zk*(1+a²), 2*(Xk+Yk*a)]          (degree 1)
// A_k = [αk+βk*a+δk*a², γk+ζk*a, εk]         (degree 2)

struct PDiffs {
    double X[4],Y[4],Z[4],S[4],P[4],Q[4];
};

static void compB(const PDiffs& d, int k, double a, double B[3]) {
    B[0]=2*(d.X[k]*a-d.Y[k]); B[1]=2*d.Z[k]*a; B[2]=-2*d.Y[k];
}
static void compD(const PDiffs& d, int k, double a, double D[2]) {
    D[0]=-2*d.Z[k]*(1+a*a); D[1]=2*(d.X[k]+d.Y[k]*a);
}
static void compA(const PDiffs& d, int k, double a, double A[3]) {
    double alpha=d.Y[k]*d.P[k]+d.Z[k]*d.Q[k];
    double beta =-d.X[k]*d.P[k]-d.Y[k]*d.S[k];
    double delta= d.X[k]*d.S[k]+d.Z[k]*d.Q[k];
    double gamma=-d.X[k]*d.Q[k]-d.Z[k]*d.S[k];
    double zeta =-d.Y[k]*d.Q[k]-d.Z[k]*d.P[k];
    double eps  = d.X[k]*d.S[k]+d.Y[k]*d.P[k];
    A[0]=alpha+beta*a+delta*a*a; A[1]=gamma+zeta*a; A[2]=eps;
}


// ── Back-substitute: solve 2x2 for (b,d) from pair (i,j) at (a,c) ────────
static bool solve_bd(const PDiffs& d, int i, int j, double a, double c,
                     double& b_out, double& d_out)
{
    double Ai[3],Bi[3],Di[2],Aj[3],Bj[3],Dj[2];
    compA(d,i,a,Ai); compB(d,i,a,Bi); compD(d,i,a,Di);
    compA(d,j,a,Aj); compB(d,j,a,Bj); compD(d,j,a,Dj);
    double Biv=peval(Bi,2,c), Div=peval(Di,1,c);
    double Bjv=peval(Bj,2,c), Djv=peval(Dj,1,c);
    double Aiv=peval(Ai,2,c), Ajv=peval(Aj,2,c);
    double det=Biv*Djv-Bjv*Div;
    if (std::abs(det)<1e-10) return false;
    b_out=(-Aiv*Djv+Ajv*Div)/det;
    d_out=(Aiv*Bjv-Ajv*Biv)/det;
    return true;
}

struct RawSol { double a, c, b, d, rsqr, eq_err; };

// ── Solve one x-parametrization ───────────────────────────────────────────
// 2D Newton on (G(a,c)=0, H(a,c)=0) from a dense 2D starting grid.
// Avoids the Sylvester resultant and its catastrophic cancellation.
static void solve_xpara(const double pts[5][3], std::vector<RawSol>& out) {
    PDiffs d;
    for (int k=1;k<=4;++k) {
        int j=k-1;
        d.X[j]=pts[k][0]-pts[0][0]; d.Y[j]=pts[k][1]-pts[0][1]; d.Z[j]=pts[k][2]-pts[0][2];
        d.S[j]=pts[k][0]+pts[0][0]; d.P[j]=pts[k][1]+pts[0][1]; d.Q[j]=pts[k][2]+pts[0][2];
    }

    // Evaluate G(a,c) = det3x3 of [rows 0,1,2], H(a,c) = det3x3 of [rows 1,2,3]
    auto eval_minor = [&](double av, double cv,
                          int r0, int r1, int r2) -> double {
        double Ar0[3],Br0[3],Dr0[2], Ar1[3],Br1[3],Dr1[2], Ar2[3],Br2[3],Dr2[2];
        compA(d,r0,av,Ar0); compB(d,r0,av,Br0); compD(d,r0,av,Dr0);
        compA(d,r1,av,Ar1); compB(d,r1,av,Br1); compD(d,r1,av,Dr1);
        compA(d,r2,av,Ar2); compB(d,r2,av,Br2); compD(d,r2,av,Dr2);
        double a0=peval(Ar0,2,cv), b0=peval(Br0,2,cv), dd0=peval(Dr0,1,cv);
        double a1=peval(Ar1,2,cv), b1=peval(Br1,2,cv), dd1=peval(Dr1,1,cv);
        double a2=peval(Ar2,2,cv), b2=peval(Br2,2,cv), dd2=peval(Dr2,1,cv);
        return a0*(b1*dd2-b2*dd1) - b0*(a1*dd2-a2*dd1) + dd0*(a1*b2-a2*b1);
    };
    auto eval_G = [&](double av, double cv){ return eval_minor(av,cv,0,1,2); };
    auto eval_H = [&](double av, double cv){ return eval_minor(av,cv,1,2,3); };

    // Back-substitute (a*,c*) → (b,d,r²), verify equations, push if accepted
    auto try_ac = [&](double a_root, double c_root) {
        int bpi=0,bpj=1; double best_det=0;
        int pairs[6][2]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
        for (auto& pr : pairs) {
            double Bi2[3],Di2[2],Bj2[3],Dj2[2];
            compB(d,pr[0],a_root,Bi2); compD(d,pr[0],a_root,Di2);
            compB(d,pr[1],a_root,Bj2); compD(d,pr[1],a_root,Dj2);
            double dt=std::abs(peval(Bi2,2,c_root)*peval(Dj2,1,c_root)
                              -peval(Bj2,2,c_root)*peval(Di2,1,c_root));
            if (dt>best_det){best_det=dt; bpi=pr[0]; bpj=pr[1];}
        }
        double b_root,d_root;
        if (!solve_bd(d,bpi,bpj,a_root,c_root,b_root,d_root)) return;

        double W=1+a_root*a_root+c_root*c_root;
        V3 off(0,b_root,d_root);
        V3 q0(pts[0][0]-off[0],pts[0][1]-off[1],pts[0][2]-off[2]);
        double qu0=q0[0]+a_root*q0[1]+c_root*q0[2];
        double rsqr=vnorm2(q0)-qu0*qu0/W;
        if (rsqr<=1e-12) return;

        double max_err=0;
        for (int k=0;k<5;++k) {
            V3 qk(pts[k][0]-off[0],pts[k][1]-off[1],pts[k][2]-off[2]);
            double quk=qk[0]+a_root*qk[1]+c_root*qk[2];
            double eq=W*(vnorm2(qk)-rsqr)-quk*quk;
            max_err=std::max(max_err,std::abs(eq));
        }
        if (max_err>1e-6*(rsqr*W)+1e-8) return;

        // Dedup within this permutation: keep best (lowest equation residual)
        for (auto& s : out) {
            if (std::abs(s.a-a_root)<1e-4*(std::abs(s.a)+std::abs(a_root)+1) &&
                std::abs(s.c-c_root)<1e-4*(std::abs(s.c)+std::abs(c_root)+1)) {
                if (max_err < s.eq_err) { s={a_root,c_root,b_root,d_root,rsqr,max_err}; }
                return;
            }
        }

        out.push_back({a_root,c_root,b_root,d_root,rsqr,max_err});
    };

    // 2D Newton from each starting point; run unconditionally, let try_ac filter
    auto run_newton = [&](double a0, double c0) {
        double a=a0, c=c0;
        double prev_res=1e300;
        int stall=0;
        for (int iter=0;iter<80;++iter) {
            double Gv=eval_G(a,c), Hv=eval_H(a,c);
            double res=std::abs(Gv)+std::abs(Hv);
            if (res>prev_res*10.) { stall++; if(stall>3) break; } else stall=0;
            prev_res=res;
            // Use finer eps as we converge
            double eps_scale=(res<1e-8)?1e-8:1e-5;
            double ea=eps_scale*(1+std::abs(a)), ec=eps_scale*(1+std::abs(c));
            double dGda=(eval_G(a+ea,c)-Gv)/ea, dGdc=(eval_G(a,c+ec)-Gv)/ec;
            double dHda=(eval_H(a+ea,c)-Hv)/ea, dHdc=(eval_H(a,c+ec)-Hv)/ec;
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

    // Starting grid: sinh-spaced in both a and c
    const int NA=40, NC=40;
    const double ASCALE=200.0, CSCALE=200.0;
    for (int ia=0;ia<NA;++ia) {
        double at=(double)ia/(NA-1)-0.5;
        double a0=ASCALE*std::sinh(5.0*at);
        for (int ic=0;ic<NC;++ic) {
            double ct=(double)ic/(NC-1)-0.5;
            double c0=CSCALE*std::sinh(5.0*ct);
            run_newton(a0,c0);
        }
    }
    // Extra starting points near axis-aligned directions (a≈0, c≈0)
    for (double av:{0.0,-0.1,0.1,-1.,1.,-5.,5.})
        for (double cv:{0.0,-0.1,0.1,-1.,1.,-5.,5.})
            run_newton(av,cv);
}

// ── Cylinder solution ──────────────────────────────────────────────────────
struct CylSol {
    V3 u_hat, offset;
    double r, phi_span, knot_err;
};

static CylSol raw_to_cyl(const RawSol& s, int perm) {
    V3 u_p(1,s.a,s.c);
    V3 o_p(0,s.b,s.d);
    double un=vnorm(u_p); V3 uh=u_p*(1.0/un);
    CylSol cs;
    // perm_map[k]: permuted coords are original[(perm_map[k][0..2])]
    // In x-para of permuted coords, axis (1,a,c) in permuted = original mapped back.
    // perm=0 (x,y,z): identity
    // perm=1 (y,z,x): perm[0]=orig[1], perm[1]=orig[2], perm[2]=orig[0]
    //   → axis (1,a,c) in (y,z,x) = original (c,1,a)
    //   un-permute: orig[0]=perm[2]=c, orig[1]=perm[0]=1, orig[2]=perm[1]=a
    //   → u_orig = (uh[2], uh[0], uh[1])
    // perm=2 (z,x,y): perm[0]=orig[2], perm[1]=orig[0], perm[2]=orig[1]
    //   → axis (1,a,c) = original (a,c,1)
    //   un-permute: orig[0]=perm[1]=a, orig[1]=perm[2]=c, orig[2]=perm[0]=1
    //   → u_orig = (uh[1], uh[2], uh[0])
    if (perm==0) {
        cs.u_hat=uh; cs.offset=o_p;
    } else if (perm==1) {
        cs.u_hat=V3(uh[2],uh[0],uh[1]);
        cs.offset=V3(o_p[2],o_p[0],o_p[1]);
    } else {
        cs.u_hat=V3(uh[1],uh[2],uh[0]);
        cs.offset=V3(o_p[1],o_p[2],o_p[0]);
    }
    cs.r=std::sqrt(s.rsqr);
    cs.phi_span=0; cs.knot_err=0;
    return cs;
}

static bool cylsols_same(const CylSol& a, const CylSol& b) {
    if (std::abs(a.r-b.r)>0.02*(a.r+b.r+1e-10)) return false;
    return std::abs(vdot(a.u_hat,b.u_hat))>0.9998;
}

static void unroll_phi(const V3 pts5[5], const CylSol& cs, double phi[5]) {
    V3 u=cs.u_hat;
    V3 tmp=(std::abs(u[0])<0.9)?V3(1,0,0):V3(0,1,0);
    V3 v=tmp-u*vdot(tmp,u); double vn=vnorm(v); v=v*(1.0/vn);
    V3 w=cross3(u,v);
    for (int k=0;k<5;++k) {
        V3 p=pts5[k]-cs.offset;
        V3 pp=p-u*vdot(p,u);
        phi[k]=std::atan2(vdot(pp,w),vdot(pp,v));
    }
    for (int k=1;k<5;++k) {
        double dv=phi[k]-phi[k-1];
        while(dv> M_PI){dv-=2*M_PI; phi[k]-=2*M_PI;}
        while(dv<-M_PI){dv+=2*M_PI; phi[k]+=2*M_PI;}
    }
}

static double knot_err_fn(const V3 pts5[5], const CylSol& cs) {
    double sum=0;
    for (int k=0;k<5;++k) {
        V3 q=pts5[k]-cs.offset;
        double qu=vdot(q,cs.u_hat);
        V3 pp=q-cs.u_hat*qu;
        double d=vnorm(pp)-cs.r; sum+=d*d;
    }
    return std::sqrt(sum/5.0);
}

static std::vector<CylSol> solve_cylinders(const V3 pts5[5]) {
    static const int perm_map[3][3]={{0,1,2},{1,2,0},{2,0,1}};
    std::vector<CylSol> all;
    for (int perm=0;perm<3;++perm) {
        double p[5][3];
        for (int i=0;i<5;++i)
            for (int d=0;d<3;++d)
                p[i][d]=pts5[i][perm_map[perm][d]];
        std::vector<RawSol> raw;
        solve_xpara(p, raw);
        for (auto& rs : raw) {
            CylSol cs=raw_to_cyl(rs,perm);
            bool dup=false;
            for (auto& ex : all) if(cylsols_same(ex,cs)){dup=true;break;}
            if (dup) continue;
            double phi[5];
            unroll_phi(pts5,cs,phi);
            cs.phi_span=std::abs(phi[4]-phi[0]);
            cs.knot_err=knot_err_fn(pts5,cs);
            all.push_back(cs);
        }
    }
    return all;
}

static double planarity_rms(const V3 pts5[5]) {
    V3 cen{};
    for (int i=0;i<5;++i) cen=cen+pts5[i];
    cen=cen*0.2;
    double M[3][3]={};
    for (int i=0;i<5;++i) {
        V3 q=pts5[i]-cen;
        for (int r=0;r<3;++r) for (int c=0;c<3;++c) M[r][c]+=q[r]*q[c];
    }
    for (int r=0;r<3;++r) for (int c=0;c<3;++c) M[r][c]/=5.0;
    fc::detail::SymEig<3> eig;
    eig.compute(M);
    return std::sqrt(std::max(0.0,eig.val[0]));
}

static void analyse_window(int win, const V3 pts5[5]) {
    double plan=planarity_rms(pts5);
    auto sols=solve_cylinders(pts5);
    int nReal=(int)sols.size();
    std::sort(sols.begin(),sols.end(),[](const CylSol&a,const CylSol&b){return a.phi_span<b.phi_span;});
    printf("win %-2d  planRMS=%.2e  nReal=%d",win,plan,nReal);
    if (nReal>0) {
        printf("  spans(deg)=[");
        for (int i=0;i<nReal;++i) printf("%.1f%s",sols[i].phi_span*180./M_PI,i+1<nReal?",":"");
        printf("]  knot=%.2e  plan/knot=%.2e  r=%.4f",
               sols[0].knot_err, plan/(sols[0].knot_err+1e-20), sols[0].r);
    } else {
        printf("  [NO REAL CYL]");
    }
    printf("\n");
}

// ── Conic evaluation ───────────────────────────────────────────────────────
static double conicV(const double co[6], double u, double vHint) {
    double A=co[0],B=co[1],C=co[2],D=co[3],E=co[4],F=co[5];
    if (std::abs(C)<1e-10) {
        double denom=B*u+E; if (std::abs(denom)<1e-10) return vHint;
        return -(A*u*u+D*u+F)/denom;
    }
    double disc=(B*u+E)*(B*u+E)-4*C*(A*u*u+D*u+F);
    if (disc<0) return vHint;
    double sq=std::sqrt(disc);
    double r1=(-(B*u+E)+sq)/(2*C), r2=(-(B*u+E)-sq)/(2*C);
    return (std::abs(r1-vHint)<std::abs(r2-vHint))?r1:r2;
}

static double pipeline_err(const V3 pts5[5], const double tvs5[5],
                            const CylSol& cs,
                            std::function<V3(double)> trueCurve)
{
    double phi5[5]; unroll_phi(pts5,cs,phi5);
    V3 u=cs.u_hat;
    V3 tmp=(std::abs(u[0])<0.9)?V3(1,0,0):V3(0,1,0);
    V3 v=tmp-u*vdot(tmp,u); v=v*(1.0/vnorm(v));
    V3 w=cross3(u,v);
    double z5[5];
    for (int k=0;k<5;++k) z5[k]=vdot(pts5[k]-cs.offset,u);
    double uvPts[5][2];
    for (int k=0;k<5;++k){uvPts[k][0]=cs.r*phi5[k]; uvPts[k][1]=z5[k];}
    double coeff[6]; fc::detail::fit_conic_2d(uvPts,coeff);
    double max_err=0;
    for (int k=0;k<4;++k) {
        for (double s:{0.25,0.5,0.75}) {
            double tMid=tvs5[k]+s*(tvs5[k+1]-tvs5[k]);
            V3 pTrue=trueCurve(tMid);
            V3 pp=pTrue-cs.offset;
            double zTrue=vdot(pp,u);
            V3 pp2=pp-u*zTrue;
            double phiTrue=std::atan2(vdot(pp2,w),vdot(pp2,v));
            double phiMid=0.5*(phi5[k]+phi5[k+1]);
            double dp=phiTrue-phiMid; dp-=2*M_PI*std::round(dp/(2*M_PI));
            phiTrue=phiMid+dp;
            double zConic=conicV(coeff,cs.r*phiTrue,0.5*(z5[k]+z5[k+1]));
            V3 pCyl=cs.offset+u*zConic+v*(cs.r*std::cos(phiTrue))+w*(cs.r*std::sin(phiTrue));
            max_err=std::max(max_err,vnorm(pCyl-pTrue));
        }
    }
    return max_err;
}

static void compare_window(int win, const V3 pts5[5], const double tvs5[5],
                            std::function<V3(double)> trueCurve)
{
    auto sols=solve_cylinders(pts5);
    int nReal=(int)sols.size();
    if (nReal<2){printf("win %-2d  nReal=%d  [can't compare]\n",win,nReal);return;}
    std::sort(sols.begin(),sols.end(),[](const CylSol&a,const CylSol&b){return a.phi_span<b.phi_span;});
    double errBest=pipeline_err(pts5,tvs5,sols[0],trueCurve);
    double errWrong=pipeline_err(pts5,tvs5,sols.back(),trueCurve);
    auto fn=[](double x)->std::string{
        if(x==0)return "0";
        int e=(int)std::floor(std::log10(std::abs(x)+1e-300));
        double m=x/std::pow(10.,e); char buf[32]; snprintf(buf,32,"%.3fe%d",m,e); return buf;
    };
    printf("win %-2d  nReal=%d  spans=[",win,nReal);
    for (int i=0;i<nReal;++i) printf("%.1f%s",sols[i].phi_span*180./M_PI,i+1<nReal?",":"");
    printf("]deg  err_correct=%s  err_wrong=%s  ratio=%s\n",
           fn(errBest).c_str(),fn(errWrong).c_str(),fn(errWrong/(errBest+1e-20)).c_str());
}

static V3 helix(double t){return V3(std::cos(2*t),std::sin(2*t),0.5*t);}
static V3 tcubic(double t){return V3(t,t*t,t*t*t);}

int main(){
    printf("=== Test A: Helix r(t)=(cos 2t, sin 2t, 0.5t), n=12 ===\n");
    printf("    exact axis=(0,0,1); WLS expects nReal=4-5, knot~1e-15\n\n");
    {
        const int n=12; std::vector<V3> ctrl(n);
        for (int i=0;i<n;++i) ctrl[i]=helix(2.0*M_PI*(double)i/(n-1));
        for (int i=0;i<n-4;++i){
            V3 pts5[5]={ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4]};
            analyse_window(i,pts5);
        }
    }
    printf("\n=== Test B: Twisted cubic r(t)=(t,t²,t³), n=12 ===\n");
    printf("    WLS expects win0 nReal=0, wins1-7 nReal=2, knot~1e-16\n\n");
    {
        const int n=12; std::vector<V3> ctrl(n); std::vector<double> tv(n);
        for (int i=0;i<n;++i){tv[i]=(double)i/(n-1); ctrl[i]=tcubic(tv[i]);}
        for (int i=0;i<n-4;++i){
            V3 pts5[5]={ctrl[i],ctrl[i+1],ctrl[i+2],ctrl[i+3],ctrl[i+4]};
            double tvs5[5]={tv[i],tv[i+1],tv[i+2],tv[i+3],tv[i+4]};
            compare_window(i,pts5,tvs5,tcubic);
        }
    }
    return 0;
}
