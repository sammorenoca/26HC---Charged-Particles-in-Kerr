#pragma once
#include <cmath>
#include <concepts>
#include <numbers>
#include <optional>

//Units G=c=M=m=1.  C++26
inline constexpr double PI = std::numbers::pi; 

struct P { double a = 0.5, E = 0.9935, L = 4.6, B = 1e-3, Q = 1e-3; };  // B=B*, Q=Q*

[[nodiscard]] constexpr double rh(const P& p) { return 1 + std::sqrt(1 - p.a * p.a); }  // horizon

// forward-mode dual number: value + gradient wrt (r, theta)
struct D {
  double v, r, t;
  constexpr D(double x = 0) : v(x), r(0), t(0) {}
  constexpr D(double x, double dr, double dt) : v(x), r(dr), t(dt) {}
};
constexpr D operator+(D x, D y) { return {x.v + y.v, x.r + y.r, x.t + y.t}; }
constexpr D operator+(D x, double y) { return {x.v + y, x.r, x.t}; }
constexpr D operator+(double x, D y) { return {x + y.v, y.r, y.t}; }
constexpr D operator-(D x, D y) { return {x.v - y.v, x.r - y.r, x.t - y.t}; }
constexpr D operator-(D x, double y) { return {x.v - y, x.r, x.t}; }
constexpr D operator-(double x, D y) { return {x - y.v, -y.r, -y.t}; }
constexpr D operator-(D x) { return {-x.v, -x.r, -x.t}; }
constexpr D operator*(D x, D y) { return {x.v * y.v, x.r * y.v + x.v * y.r, x.t * y.v + x.v * y.t}; }
constexpr D operator*(D x, double y) { return {x.v * y, x.r * y, x.t * y}; }
constexpr D operator*(double x, D y) { return {x * y.v, x * y.r, x * y.t}; }
constexpr D operator/(D x, D y) {
  const double i = 1 / y.v, q = x.v * i;
  return {q, (x.r - q * y.r) * i, (x.t - q * y.t) * i};
}
constexpr D operator/(double x, D y) {
  const double i = 1 / y.v, q = x * i;
  return {q, -q * y.r * i, -q * y.t * i};
}
inline D sin(D x) { return {std::sin(x.v), std::cos(x.v) * x.r, std::cos(x.v) * x.t}; }
inline D cos(D x) { return {std::cos(x.v), -std::sin(x.v) * x.r, -std::sin(x.v) * x.t}; }

template <class T> concept Scalar = std::same_as<T, double> || std::same_as<T, D>;

// Eq.(6): H1(r,theta) with metric Eq.(1), potentials Eqs.(3)-(4), f1..f3 Eqs.(7)-(9)
template <Scalar T> [[gnu::always_inline]] inline T H1(T r, T th, const P& p) {
  using std::sin; using std::cos;
  const T s = sin(th), c = cos(th), s2 = s * s;
  const T Sg = r * r + p.a * p.a * c * c;              // Sigma
  const T w = 2.0 * r * (1.0 / Sg);
  const T gtt = -(1.0 - w);
  const T gtp = -w * p.a * s2;
  const T gpp = (r * r + p.a * p.a + w * p.a * p.a * s2) * s2;
  const T k = p.a * p.B - 0.5 * p.Q;
  const T At = k * gtt + 0.5 * p.B * gtp;              // q*A_t
  const T Ap = k * gtp + 0.5 * p.B * gpp;              // q*A_phi
  const T d = gtt * gpp - gtp * gtp;
  const T id = 1.0 / d;
  const T f1 = gpp * id, f2 = gtp * id, f3 = gtt * id; // Eqs.(7)-(9)
  const T X = p.E + At, Y = p.L - Ap;
  const T u = f1 * X + f2 * Y, v = f2 * X + f3 * Y;
  return 0.5 * gtt * u * u + 0.5 * gpp * v * v - gtp * u * v;
}

// Eq.(13): effective potential. H1 is quadratic in E -> A,B,C from 3 samples.
[[nodiscard]] inline double veff(double r, const P& p) {
  P q = p;
  q.E = 0;  const double g = H1(r, PI / 2, q);
  q.E = 1;  const double hp = H1(r, PI / 2, q);
  q.E = -1; const double hm = H1(r, PI / 2, q);
  const double A = g - 0.5 * (hp + hm), B = 0.5 * (hp - hm), C = g;
  return B / (2 * A) + std::sqrt(B * B + 4 * A * C + 2 * A) / (2 * A);
}

// Eqs.(12)-(13): H = -1/2 at theta=pi/2, pr=0 -> p_theta (nullopt if the orbit is forbidden)
[[nodiscard]] inline std::optional<double> init_pt(double r, const P& p) {
  const double v = r * r * (-1 - 2 * H1(r, PI / 2, p));
  if (!(v > 0)) return std::nullopt;
  return std::sqrt(v);
}
