#include "inmost.h"
#include "solver.cpp"
#include <iostream>
#include <functional>
#include <cmath>
#include <chrono>

using namespace INMOST;

double f(double x, double y)
{
    return 26.0 * sin(5 * x) * cos(y);
}

double u(double x, double y)
{
    return sin(5 * x) * cos(y);
}

double norm_L2_01(INMOST::Sparse::Vector& x, std::function<double(double, double)> u, size_t n)
{
    double h = 1.0 / n;
    double norm = 0.0;
    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 1; j < n; ++j) {
            double diff = u(i * h, j * h) - x[(i - 1) * (n - 1) + j - 1];
            norm += diff * diff;
        }
    }
    return sqrt(norm);
}

double norm_C_01(INMOST::Sparse::Vector& x, std::function<double(double, double)> u, size_t n)
{
    double h = 1.0 / n;
    double norm = 0.0;
    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 1; j < n; ++j) {
            norm = std::max(norm, abs(u(i * h, j * h) - x[(i - 1) * (n - 1) + j - 1]));
        }
    }
    return sqrt(norm);
}

int main(int argc, char *argv[]) {
    size_t max_n = 600;
	
    double dx = 1.0, dy = 1.0;
    for (size_t n = 2; n < max_n; n *= 2) {
        std::cout << n << ' ';
        auto start = std::chrono::steady_clock::now();
        Sparse::Vector x = solve_DE(f,
                    [] (double x) { return u(x, 0.0); },
                    [] (double x) { return u(x, 1.0); },
                    [] (double y) { return u(0.0, y); },
                    [] (double y) { return u(1.0, y); }, 
		    dx, dy, n);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> time = end - start;
        std::cout << time.count() << ' ' << norm_L2_01(x, u, n) << ' ' << norm_C_01(x, u, n) << std::endl;
    }
    return 0;
}
