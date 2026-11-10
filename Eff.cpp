#include "kerr.hpp"
#include <array>
#include <fstream>
#include <format>
#include <numeric>
#include <print>
#include <span>
#include <string>
#include <string_view>
#include <vector>

struct Isco { double R, L, E; };

double d1(const P& p, double r) { 
  constexpr double h = 1e-4; 
  return (veff(r + h, p) - veff(r - h, p)) / (2 * h); 
}
double d2(const P& p, double r) { 
  constexpr double h = 2e-3;
  return (veff(r + h, p) - 2 * veff(r, p) + veff(r - h, p)) / (h * h); 
}

// stable circular orbits: dV/dr=0 with V''>0, Eqs.(14)-(15)
[[nodiscard]] std::vector<double> stab(const P& p, double hi = 80, double dr = 0.02) {
  std::vector<double> o;
  double r = rh(p) + 0.2, f = d1(p, r);
  
  for (; r < hi; r += dr) {
    const double g = d1(p, r + dr);
    
    if (f < 0 && g > 0) {
      double a = r, b = r + dr;
      for (int i = 0; i < 60; ++i) { 
        const double m = std::midpoint(a, b); 
        (d1(p, m) < 0 ? a : b) = m; 
      }
      o.push_back(std::midpoint(a, b));
    }
    f = g;
  }
  return o;
}

// ISCO: V'=V''=0 (equality in Eq.15), Newton on (R,L); E=V(R)
[[nodiscard]] std::optional<Isco> isco(P p, double R, double L) {
  const auto F = [&](double r, double l) { P q = p; q.L = l; return std::pair{d1(q, r), d2(q, r)}; };
  for (int i = 0; i < 80; ++i) {
    constexpr double hr = 1e-3, hl = 1e-4;
    const auto [a1, a2] = F(R + hr, L);
    const auto [b1, b2] = F(R - hr, L);
    const auto [c1, c2] = F(R, L + hl);
    const auto [e1, e2] = F(R, L - hl);
    const auto [f1, f2] = F(R, L);
    const double j11 = (a1 - b1) / (2 * hr), j21 = (a2 - b2) / (2 * hr);
    const double j12 = (c1 - e1) / (2 * hl), j22 = (c2 - e2) / (2 * hl);
    const double det = j11 * j22 - j12 * j21;
    const double dR = (f1 * j22 - f2 * j12) / det, dL = (j11 * f2 - j21 * f1) / det;
    R -= dR; L -= dL;
    if (std::fabs(dR) + std::fabs(dL) < 1e-10) break;
  }
  p.L = L;
  const double E = veff(R, p);
  if (!std::isfinite(E)) return std::nullopt;
  return Isco{R, L, E};
}

// analytic Kerr ISCO (Bardeen), s=+1 direct, -1 retrograde: check for B=Q=0
[[nodiscard]] double kerr_isco(double a, int s) {
  const double z1 = 1 + std::cbrt(1 - a * a) * (std::cbrt(1 + a) + std::cbrt(1 - a));
  const double z2 = std::sqrt(3 * a * a + z1 * z1);
  return 3 + z2 - s * std::sqrt((3 - z1) * (3 + z1 + 2 * z2));
}

[[nodiscard]] std::string lab(std::string_view k, double x) { return std::format("{}={:g}", k, x); }

void dumpV(const std::string& fn, std::span<const P> ps, double r1, std::span<const std::string> lb) {
  std::ofstream f(fn);
  std::print(f, "r");
  for (const auto& l : lb) std::print(f, ",{}", l);
  std::println(f);
  for (double r = rh(ps[0]) + 0.05; r <= r1; r += 0.05) {
    std::print(f, "{:.3f}", r);
    for (const auto& p : ps) std::print(f, ",{:.8f}", veff(r, p));
    std::println(f);
  }
}

void figs() {  // Figs.1-4 data
  constexpr double L0 = 2 * 1.7320508075688772;  // 2*sqrt(3)
  std::vector<P> v;
  std::vector<std::string> l;
  const auto add = [&](const P& p, std::string s) { v.push_back(p); l.push_back(std::move(s)); };
  const auto flush = [&](const std::string& fn, double r1) { dumpV(fn, v, r1, l); v.clear(); l.clear(); };

  for (int i = 1; i <= 5; ++i) { const double b = 1e-3 * i; add({.a = 0.1, .L = L0, .B = b, .Q = 2e-4}, lab("B", b)); }
  flush("fig1_veff_B.csv", 20);
  for (const double q : {0.0, 1e-4, 4e-4}) add({.a = 0.1, .L = L0, .B = 1e-3, .Q = q}, lab("Q", q));
  flush("fig2_veff_Q.csv", 10);
  for (const double a : {0.02, 0.04, 0.06, 0.08, 0.1}) add({.a = a, .L = L0, .B = 1e-3, .Q = 4e-3}, lab("a", a));
  flush("fig3_veff_a.csv", 15);
  const std::array<std::pair<double, std::array<double, 4>>, 3> f4 = {{
      {0.0, {0, 3.4641, 4, 4.5}}, {0.4, {0, 3.0340, 4, 4.5}}, {0.6, {0, 2.7559, 2.8, 2.9}}}};
  for (const auto& [a, Ls] : f4) {
    for (const double L : Ls) add({.a = a, .L = L, .B = 0, .Q = 0}, lab("L", L));
    flush(std::format("fig4_veff_L_a{:.1f}.csv", a), 50);
  }
}

void tab12() {
  constexpr double L0 = 2 * 1.7320508075688772;
  constexpr std::array<double, 8> r1 = {7.4595, 7.5165, 7.5693, 7.6177, 7.6618, 7.7013, 7.7365, 7.7673};
  std::println("\n== Table 1 (a=0.1, L=2sqrt3, Q*=2e-4): R of stable orbit ==\n  B*      R(here)   R(paper)");
  for (int i = 0; i < 8; ++i) {
    const P p{.a = 0.1, .L = L0, .B = 1e-3 * (i + 1), .Q = 2e-4};
    const auto s = stab(p);
    std::println("  {:.3f}   {:.4f}    {:.4f}", p.B, s.empty() ? NAN : s[0], r1[i]);
  }
  constexpr double r2[3][2] = {{7.4559, 7.5129}, {7.4577, 7.5147}, {7.4631, 7.5201}};
  constexpr double qs[3] = {0, 1e-4, 4e-4};
  std::println("\n== Table 2 (a=0.1): R vs Q* ==\n  Q*       B*=0.001 (paper)   B*=0.002 (paper)");
  for (int i = 0; i < 3; ++i) {
    std::print("  {:.4f} ", qs[i]);
    for (int j = 0; j < 2; ++j) {
      const P p{.a = 0.1, .L = L0, .B = 1e-3 * (j + 1), .Q = qs[i]};
      std::print("  {:.4f} ({:.4f})   ", stab(p)[0], r2[i][j]);
    }
    std::println();
  }
}

void tab3() {
  struct C { double a, L, R; };
  constexpr std::array<C, 8> c = {{{0, 4, 12}, {0, 4.5, 16.5876}, {0.2, 4, 12.5329}, {0.2, 4.5, 16.9761},
                                   {0.4, 4, 12.9954}, {0.4, 4.5, 17.3336}, {0.6, 2.8, 4.5721}, {0.6, 2.9, 5.4096}}};
  std::println("\n== Table 3 (B*=Q*=0): stable orbit radius ==\n  a    L     R(here)   R(paper)");
  for (const auto& x : c) {
    const P p{.a = x.a, .L = x.L, .B = 0, .Q = 0};
    std::println("  {:.1f}  {:.1f}   {:.4f}   {:.4f}", x.a, x.L, stab(p)[0], x.R);
  }
  constexpr std::array<C, 4> ic = {{{0, 3.4641, 6}, {0.2, 3.2640, 5.3294}, {0.4, 3.0340, 4.6143}, {0.6, 2.7559, 3.8290}}};
  std::println("  ISCO: a   R(here)  L(here)   R(paper) L(paper)  R(Bardeen)");
  for (const auto& x : ic) {
    const auto r = isco({.a = x.a, .B = 0, .Q = 0}, x.R, x.L);
    if (r) std::println("        {:.1f}  {:.4f}  {:.4f}    {:.4f}  {:.4f}   {:.4f}", x.a, r->R, r->L, x.R, x.L, kerr_isco(x.a, 1));
  }
}

// Tables 4-6: ISCO (R,L,E) direct / retrograde
void tab456() {
  struct C { std::string_view t; double a, B, Q; std::array<std::array<double, 3>, 2> ref; };
  const std::array<C, 9> c = {{
      {"T4 B*=0.02 Q*=0.008", 0.2, 0.02, 0.008, {{{5.2406, 3.1794, 0.9078}, {6.4737, -3.8201, 0.9879}}}},
      {"T4 B*=0.04 Q*=0.008", 0.2, 0.04, 0.008, {{{5.0411, 3.1433, 0.8885}, {6.1176, -4.0808, 1.0406}}}},
      {"T4 B*=0.06 Q*=0.008", 0.2, 0.06, 0.008, {{{4.7844, 3.1311, 0.8726}, {5.8262, -4.3879, 1.1017}}}},
      {"T5 B*=0.01 Q*=0    ", 0.2, 0.01, 0.0, {{{5.3276, 3.2224, 0.9226}, {6.5990, -3.7282, 0.9688}}}},
      {"T5 B*=0.01 Q*=0.001", 0.2, 0.01, 0.001, {{{5.3275, 3.2209, 0.9222}, {6.6193, -3.7265, 0.9683}}}},
      {"T5 B*=0.01 Q*=0.002", 0.2, 0.01, 0.002, {{{5.3274, 3.2194, 0.9217}, {6.6043, -3.7247, 0.9679}}}},
      {"T6 B*=0    Q*=0    ", 0.2, 0.0, 0.0, {{{5.3641, 3.2641, 0.9354}, {6.6713, -3.6434, 0.9485}}}},
      {"T6 B*=0.02 Q*=0    ", 0.2, 0.02, 0.0, {{{5.2534, 3.1912, 0.9112}, {6.4653, -3.8343, 0.9917}}}},
      {"T6 B*=0.04 Q*=0    ", 0.2, 0.04, 0.0, {{{5.0411, 3.1546, 0.8918}, {6.1218, -4.0954, 1.0445}}}}}};
  std::println("\n== Tables 4-6: ISCO  (here | paper)  R, L, E ==");
  for (const auto& x : c)
    for (int k = 0; k < 2; ++k) {
      const auto& rf = x.ref[k];
      const auto r = isco({.a = x.a, .B = x.B, .Q = x.Q}, rf[0], rf[1]);
      if (r) std::println("  {} {}  R {:.4f}|{:.4f}  L {:.4f}|{:.4f}  E {:.4f}|{:.4f}", x.t, k ? "RO" : "DO", r->R, rf[0], r->L, rf[1], r->E, rf[2]);
    }
}

// NEW 1: B=Q=0 against analytic Kerr ISCO
void new_bardeen() {
  std::println("\n== NEW 1: ISCO vs analytic Kerr (B*=Q*=0) ==\n  a     s   R(num)    R(Bardeen)  |diff|");
  for (const double a : {0.0, 0.2, 0.5, 0.9})
    for (const int s : {1, -1}) {
      const double ref = kerr_isco(a, s);
      const auto r = isco({.a = a, .B = 0, .Q = 0}, ref + 0.05, s * 3.0);
      if (r) std::println("  {:.1f}  {:+d}  {:.8f}  {:.8f}  {:.1e}", a, s, r->R, ref, std::fabs(r->R - ref));
    }
}

// NEW 2: ISCO shift and binding energy vs B*, with Wald charge Q*=2aB* and Q*=0
void new_wald() {
  std::println("\n== NEW 2: a=0.2 ISCO vs B*; Q*=0 vs Wald Q*=2aB* ==\n  B*     Q*      sense   R        L        E        1-E");
  std::ofstream f("isco_wald.csv");
  std::println(f, "B,Q,sense,R,L,E");
  const std::array<std::array<double, 2>, 2> g0 = {{{5.36, 3.26}, {6.67, -3.64}}};  // guesses at B=0
  auto gr = std::array{g0, g0};                                                    // continuation guesses [wald][sense]
  for (const double B : {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06})
    for (int w = 0; w < 2; ++w) {
      if (B == 0 && w) continue;
      const double a = 0.2, Q = w ? 2 * a * B : 0;
      for (int k = 0; k < 2; ++k) {
        const auto& g = (B == 0) ? g0[k] : gr[w][k];
        const auto r = isco({.a = a, .B = B, .Q = Q}, g[0], g[1]);
        if (!r) continue;
        gr[w][k] = {r->R, r->L};
        std::println("  {:.2f}  {:.4f}  {}    {:.4f}  {:+.4f}  {:.4f}  {:+.4f}", B, Q, k ? "RO" : "DO", r->R, r->L, r->E, 1 - r->E);
        std::println(f, "{:g},{:g},{},{:.6f},{:.6f},{:.6f}", B, Q, k, r->R, r->L, r->E);
      }
    }
}

int main() {
  figs(); tab12(); tab3(); tab456(); new_bardeen(); new_wald();
}
