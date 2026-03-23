// diag_parabola.cpp — diagnose which parabola windows fail ConicWindow::valid()
//
// Compile: g++ -std=c++17 -O2 -I. -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 -o diag_parabola diag_parabola.cpp
// Run:     ./diag_parabola

#define FC_DEBUG_PARABOLA 1
#include "conicblend.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

using namespace fc;
using V2 = VecN<2>;

static constexpr double PI = 3.14159265358979323846;

static std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i) v[i] = a + (b - a) * i / (n - 1);
    return v;
}

int main() {
    int n = 14;
    auto ts = linspace(-2.0, 2.0, n);
    std::vector<V2> ctrl(n);
    for (int i = 0; i < n; ++i) ctrl[i] = V2(ts[i], ts[i]*ts[i]);

    printf("Parabola y=x², n=%d, %d windows\n", n, n-4);
    printf("%-4s %-8s %-8s  %s\n", "win", "t_lo", "t_hi", "valid");
    printf("%s\n", std::string(40, '-').c_str());

    int pass = 0, fail = 0;
    for (int i = 0; i < n - 4; ++i) {
        ConicWindow<2> w(ctrl[i], ctrl[i+1], ctrl[i+2], ctrl[i+3], ctrl[i+4],
                         ts[i], ts[i+1], ts[i+2], ts[i+3], ts[i+4]);
        bool v = w.valid();
        printf("%-4d %-8.4f %-8.4f  %s\n", i, ts[i], ts[i+4], v ? "PASS" : "FAIL");
        if (v) ++pass; else ++fail;
    }
    printf("\n%d/%d pass\n", pass, pass+fail);
    return 0;
}
