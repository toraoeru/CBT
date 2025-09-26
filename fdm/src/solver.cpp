#include <vector>
#include <functional>
#include <iostream>
#include <cmath>
#include "inmost.h"

using namespace INMOST;

Sparse::Vector
solve_DE(std::function<double(double, double)> f,
            std::function<double(double)> g_bottom,
            std::function<double(double)> g_top,
            std::function<double(double)> g_left,
            std::function<double(double)> g_right, 
	    double dx, double dy, size_t n)
{
    double h = 1.0 / n;
    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector x;

    A.SetInterval(0, (n - 1) * (n - 1));
    b.SetInterval(0, (n - 1) * (n - 1));
    x.SetInterval(0, (n - 1) * (n - 1));

    std::function<void(size_t, size_t, size_t, bool)> set_coeff =
        [&] (size_t row, size_t i_set, size_t j_set, bool mode) {
            if (i_set == 0) {
                b[row] += g_bottom(j_set * h) / h / h;
            } else if (i_set == n) {
                b[row] += g_top(j_set * h) / h / h;
            } else if (j_set == 0) {
                b[row] += g_left(i_set * h) / h / h;
            } else if (j_set == n) {
                b[row] += g_right(i_set * h) / h / h;
            } else {
                A[row][(i_set - 1) * (n - 1) + j_set - 1] = mode ? -dx : -dy;
            }
        };

    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 1; j < n; ++j) {
            size_t idx = (i - 1) * (n - 1) + j - 1;
            A[idx][idx] = 2.0 * (dx + dy);
            b[idx] = f(i * h, j * h);
            set_coeff(idx, i - 1, j, true);
            set_coeff(idx, i + 1, j, true);
            set_coeff(idx, i, j - 1, false);
            set_coeff(idx, i, j + 1, false);
        }
    }

    Solver S(Solver::INNER_DDPQILUC);
    S.SetParameter("absolute_tolerance", "1e-12");
    S.SetParameter("relative_tolerance", "1e-12");
    S.SetParameter("drop_tolerance", "0.005");
    S.SetMatrix(A);
    S.Solve(b, x);
    std::cout << S.Iterations() << ' ' << S.IterationsTime() << ' ';

    for (size_t i = 0; i < (n - 1) * (n - 1); ++i) {
        x[i] *= h * h;
    }
    return x;
}
