#ifndef KERR_HPP
#define KERR_HPP
 
// Unidades geometrizadas: c = G = M = m = 1.  Q* = qQ, B* = qB.
// Referencia: Sun, Wu et al., arXiv:2111.04900v1

#include <cmath>
 
namespace kerr {
 
struct Params {
    double a;   
    double Qs;  
    double Bs;  
    double L;   
};
 
struct Metric {
    double gtt;
    double gtp; 
    double gpp;  
};
 
inline Metric metric(double r, double th, double a) {

    const double s2  = std::sin(th) * std::sin(th);
    const double Sig = r * r + a * a * std::cos(th) * std::cos(th);
    const double rho2 = r * r + a * a;

    Metric g;
    g.gtt = -(1.0 - 2.0 * r / Sig);
    g.gtp = -(2.0 * r * a * s2) / Sig;
    g.gpp =  (rho2 + 2.0 * r * a * a * s2 / Sig) * s2;
    
    return g;
}
 
struct FFunc {
    double f1, f2, f3;
};
 
inline FFunc ffunc(const Metric& g) {

    const double det = g.gtt * g.gpp - g.gtp * g.gtp;

    FFunc f;
    f.f1 = g.gpp / det;
    f.f2 = g.gtp / det;
    f.f3 = g.gtt / det;
    return f;
}
 
struct Coef {
    double A, B, C;
};
 
inline Coef coef(double r, double th, const Params& p) {

    const Metric g = metric(r, th, p.a);
    const FFunc  f = ffunc(g);

    const double f1 = f.f1, f2 = f.f2, f3 = f.f3;
    const double gtt = g.gtt, gtp = g.gtp, gpp = g.gpp;
    const double L = p.L;
 
  
    const double At = (p.a * p.Bs - 0.5 * p.Qs) * gtt + 0.5 * p.Bs * gtp;    // A_t   = (a B* - Q*/2) g_tt + (B*/2) g_tphi
    const double Ap = (p.a * p.Bs - 0.5 * p.Qs) * gtp + 0.5 * p.Bs * gpp;    // A_phi = (a B* - Q*/2) g_tphi + (B*/2) g_phiphi
 
    const double A = -0.5 * (f1 * f1 * gtt + f2 * f2 * gpp - 2.0 * f1 * f2 * gtp);
 
    const double B1 = gtt * At * f1 * f1 + gtt * f1 * f2 * L - gtt * f1 * f2 * Ap;
    const double B2 = gpp * At * f2 * f2 + gpp * f2 * f3 * L - gpp * f2 * f3 * Ap;
    const double B3 = -2.0 * gtp * At * f1 * f2 - gtp * f1 * f3 * L + gtp * Ap * f1 * f3
                      - gtp * f2 * f2 * L + gtp * Ap * f2 * f2;
    const double B = B1 + B2 + B3;
 
    const double C1 = f1 * f1 * At * At + 2.0 * f1 * f2 * L * At - 2.0 * f1 * f2 * At * Ap
                      + f2 * f2 * (Ap * Ap - 2.0 * Ap * L + L * L);
    const double C2 = f2 * f2 * At * At + 2.0 * f2 * f3 * L * At - 2.0 * f2 * f3 * At * Ap
                      + f3 * f3 * (Ap * Ap - 2.0 * Ap * L + L * L);
    const double C3 = Ap * Ap * f2 * f3 - Ap * At * f1 * f3 - Ap * At * f2 * f2
                      + At * At * f1 * f2 - 2.0 * Ap * f2 * f3 * L + At *
                      + At * f2 * f2 * L + f2 * f3 * L * L;

    const double C = 0.5 * gtt * C1 + 0.5 * gpp * C2 - gtp * C3;
 
    return Coef{A, B, C};
}
 
inline double Veff(double r, double th, const Params& p) {
    const Coef c = coef(r, th, p);
    const double disc = c.B * c.B + 4.0 * c.A * c.C + 2.0 * c.A;

    return c.B / (2.0 * c.A) + std::sqrt(disc) / (2.0 * c.A); // V = E = B/(2A) + sqrt(B^2 + 4AC + 2A)/(2A)
}
 
inline double rh(double a) { 
    return 1.0 + std::sqrt(1.0 - a * a); 
} 
 
} 
 
#endif 