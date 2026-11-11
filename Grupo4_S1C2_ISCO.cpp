#include "Grupo4_S1C2_ISCO.hpp"
#include <Eigen/Dense>
#include <numbers>
#include <cmath>
#include <cstdio>

typedef Eigen::Vector2d Vec;
double h = 1e-6;

double V(double r, kerr::Params p) { 
    return kerr::Veff(r, std::numbers::pi / 2.0, p);          // Eq. (13)                        
}   

double dV(double r, kerr::Params p) {   // dV/dr
    return (-V(r + 2*h, p) + 8*V(r + h, p) - 8*V(r - h, p) + V(r - 2*h, p)) / (12*h);
}

double d2V(double r, kerr::Params p) {   // d2V/dr2
    double h = 1e-6;
    return (-V(r + 2*h, p) + 16*V(r + h, p) - 30*V(r, p) + 16*V(r - h, p) - V(r - 2*h, p)) / (12*h*h);
}

double stable(kerr::Params p) {

    double R = 0.0;
    double r0 = kerr::rh(p.a) + 0.06;
    double f0 = dV(r0, p);

    for (double r = r0; r < 60; r += 0.005) {
        double f1 = dV(r + 0.005, p);

        if (f0 < 0 && f1 >= 0){
            R = r;
        }

        f0 = f1;
    }

    for (int i = 0; i < 5; i++){
        R -= dV(R, p) / d2V(R, p);
    }

    return R;
}

// F(r, L) = (dV/dr, d2V/dr2) = 0  (Eqs. 14-15 with equality)
Vec F(Vec x, kerr::Params p) {
    p.L = x[1];
    return Vec(dV(x[0], p), d2V(x[0], p));
}

Vec ISCO(kerr::Params p, Vec x) {   // x = (R, L),  Multidimensional Newton Raphson
    double d = 1e-3;
    for (int i = 0; i < 30; i++) {
        Eigen::Matrix2d J;
        J.col(0) = (F(x + Vec(d, 0), p) - F(x - Vec(d, 0), p)) / (2*d);
        J.col(1) = (F(x + Vec(0, d), p) - F(x - Vec(0, d), p)) / (2*d);
        x -= J.colPivHouseholderQr().solve(F(x, p));
    }
    return x;
}

void rows(FILE* t, double x, kerr::Params p, Vec g[2]) {   
    for (int s = 0; s < 2; s++) {
        g[s] = ISCO(p, g[s]);
        p.L = g[s][1];
        fprintf(t, "%g,%s,%.4f,%.4f,%.5f\n", x, s ? "RO" : "DO", g[s][0], g[s][1], V(g[s][0], p));
    }
}

void curve(FILE* f, const char* id, double val, kerr::Params p, double r0, double r1) {
    r0 = std::fmax(r0, kerr::rh(p.a) + 0.01);
    for (int i = 0; i <= 1500; i++) {
        double r = r0 + (r1 - r0) * i / 1500;
        fprintf(f, "%s,%.12g,%.12g,%.12g\n", id, val, r, V(r, p));
    }
}

int main() {
    FILE* t = fopen("tablas.txt", "w");
    double L0 = 2 * std::sqrt(3.0);

    fprintf(t, "# Table 1: a=0.1, Q=2e-4, L=2sqrt(3)\nB,R\n");
    for (int k = 1; k <= 8; k++){
        fprintf(t, "%g,%.4f\n", k * 1e-3, stable({0.1, 2e-4, k * 1e-3, L0}));
    }

    fprintf(t, "\n# Table 2: a=0.1, L=2sqrt(3)\nQ,B,R\n");
    for (double Q : {0.0, 1e-4, 4e-4})
        for (double B : {1e-3, 2e-3}){
            fprintf(t, "%g,%g,%.4f\n", Q, B, stable({0.1, Q, B, L0}));
        }

    fprintf(t, "\n# Table 3: Q=B=0\na,L,R,type\n");

    double a3[4] = {0, 0.2, 0.4, 0.6};
    double L3[4][2] = {{4, 4.5}, {4, 4.5}, {4, 4.5}, {2.8, 2.9}};
    double Lc[4];

    Vec x(6, 3.46);
    for (int i = 0; i < 4; i++) {
        x = ISCO({a3[i], 0, 0, 0}, x);
        Lc[i] = x[1];

        fprintf(t, "%g,%.4f,%.4f,ISCO\n", a3[i], x[1], x[0]);

        for (double L : L3[i]){
            fprintf(t, "%g,%.4f,%.4f,stable\n", a3[i], L, stable({a3[i], 0, 0, L}));
        }
    }

    Vec g[2] = {Vec(5.3294433, 3.2640286), Vec(6.63904, -3.6433589)};
    fprintf(t, "\n# Table 4: a=0.2, Q=0.008\nB,orbit,R,L,E\n");
    for (double B : {0.02, 0.04, 0.06}){
         rows(t, B, {0.2, 0.008, B, 0}, g);

    }

    g[0] = Vec(5.3294433, 3.2640286); 
    g[1] = Vec(6.63904, -3.6433589);

    fprintf(t, "\n# Table 5: a=0.2, B=0.01\nQ,orbit,R,L,E\n");
    for (double Q : {0.0, 0.001, 0.002}){
        rows(t, Q, {0.2, Q, 0.01, 0}, g);
    }

    g[0] = Vec(5.3294433, 3.2640286); 
    g[1] = Vec(6.63904, -3.6433589);

    fprintf(t, "\n# Table 6: a=0.2, Q=0\nB,orbit,R,L,E\n");
    for (double B : {0.0, 0.02, 0.04}){
        rows(t, B, {0.2, 0, B, 0}, g);
    }

    fclose(t);

    FILE* c = fopen("curvas.csv", "w");
    fprintf(c, "panel,valor,r,V\n");
    for (int k = 1; k <= 5; k++) curve(c, "fig1", k * 1e-3, {0.1, 2e-4, k * 1e-3, L0}, 0, 20);
    for (double Q : {0.0, 1e-4, 4e-4}) curve(c, "fig2", Q, {0.1, Q, 1e-3, L0}, 3, 10);
    for (double a : {0.02, 0.04, 0.06, 0.08, 0.10}) curve(c, "fig3", a, {a, 0.004, 1e-3, L0}, 0, 15);
    for (double L : {0.0, Lc[0], 4.0, 4.5}) curve(c, "fig4a", L, {0, 0, 0, L}, 0, 50);
    for (double L : {0.0, Lc[2], 4.0, 4.5}) curve(c, "fig4b", L, {0.4, 0, 0, L}, 0, 50);
    for (double L : {0.0, Lc[3], 2.8, 2.9}) curve(c, "fig4c", L, {0.6, 0, 0, L}, 0, 35);
    fclose(c);
}