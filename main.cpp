#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric>

std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] + b[i];
    }
    return result;
}

double F(double dy, double x, double y) {
    return 2 * x + std::sin(dy) - std::cos(y);
}

double Fy(double dy, double x, double y) {
    return std::sin(y);
}

double Fdy(double dy, double x, double y) {
    return std::cos(dy);
}

void tridag(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& r, std::vector<double>& u) {
    int n = a.size();
    std::vector<double> gam(n);
    if (b[0] == 0.0) throw std::runtime_error("Error: Division by zero in tridag.");
    double bet = b[0];
    u[0] = r[0] / bet;

    for (int j = 1; j < n; j++) {
        gam[j] = c[j - 1] / bet;
        bet = b[j] - a[j] * gam[j];
        if (bet == 0.0) throw std::runtime_error("Error: Division by zero during tridag computation.");
        u[j] = (r[j] - a[j] * u[j - 1]) / bet;
    }
    for (int j = n - 2; j >= 0; j--) {
        u[j] -= gam[j + 1] * u[j + 1];
    }
}

void JCalc(double alpha, double beta, const std::vector<double>& x, const std::vector<double>& y, int n, double h, bool firstRun,
           std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& r) {
    a.assign(n, 0.0);
    b.assign(n, 0.0);
    c.assign(n, 0.0);
    r.assign(n, 0.0);
    
    if (firstRun) {
        b[0] = 2 + std::pow(h, 2) * Fy((beta - alpha) / (2 * h), x[0], y[0]);
        r[0] = alpha - 2 * y[0] + beta - std::pow(h, 2) * F((beta - alpha) / (2 * h), x[0], y[0]);
    } else {
        b[0] = 2 + std::pow(h, 2) * Fy((y[1] - alpha) / (2 * h), x[0], y[0]);
        r[0] = alpha - 2 * y[0] + y[1] - std::pow(h, 2) * F((y[1] - alpha) / (2 * h), x[0], y[0]);
        
        for (int j = 1; j < n - 1; j++) {
            double dy = (y[j + 1] - y[j - 1]) / (2 * h);
            a[j] = -1 - (h / 2) * Fdy(dy, x[j], y[j]);
            b[j] = 2 + std::pow(h, 2) * Fy(dy, x[j], y[j]);
            c[j] = -1 + (h / 2) * Fdy(dy, x[j], y[j]);
            r[j] = y[j - 1] - 2 * y[j] + y[j + 1] - std::pow(h, 2) * F(dy, x[j], y[j]);
        }
    }
    a[n - 1] = -1 - (h / 2) * Fdy((beta - y[n - 2]) / (2 * h), x[n - 1], y[n - 1]);
    b[n - 1] = 2 + std::pow(h, 2) * Fy((beta - y[n - 2]) / (2 * h), x[n - 1], y[n - 1]);
    r[n - 1] = y[n - 2] - 2 * y[n - 1] + beta - std::pow(h, 2) * F((beta - y[n - 2]) / (2 * h), x[n - 1], y[n - 1]);
}

std::vector<double> yUpdater(double a, double b, double alpha, double beta, std::vector<double>& y, const std::vector<double>& x, int n, double h) {
    std::vector<double> yPrev = y;
    y.resize(n);
    if (y.size() == 1) {
        return yPrev;
    }
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            y[i] = (yPrev[0] + alpha) / 2;
        } else if (i == n - 1) {
            y[i] = (beta + yPrev[i / 2 - 1]) / 2;
        } else if (i % 2 == 0) {
            y[i] = (yPrev[i / 2] + yPrev[i / 2]) / 2;
        } else {
            y[i] = yPrev[i / 2];
        }
    }
    return y;
}

std::vector<double> xUpdater(double a, double b, double alpha, double beta, const std::vector<double>& y, std::vector<double>& x, int n, double h) {
    x.resize(n);
    for (int i = 0; i < n; i++) {
        x[i] = a + (i + 1) * h;
    }
    return x;
}

void finiteDiff(double alpha, double beta, double x0, double xEnd, double accuracy) {
    std::vector<double> prev_y(1, alpha), prev_prev_y(1, alpha);
    double targetIndex = 0, prev_targetIndex = 0, prev_prev_targetIndex = 0;

    std::cout << std::setw(2) << "i"; 
    std::cout << std::setw(15) << "A(h_i)";
    std::cout << std::setw(20) << "A(h_(i-1))-A(h_i)"; 
    std::cout << std::setw(15) << "alpha_k";
    std::cout << std::setw(15) << "rich-error"; 
    std::cout << std::setw(15) << "order";
    std::cout << std::setw(10) << "n";
    std::cout << std::endl;

    std::vector<double> y(1, 0.5 * (alpha + beta));
    std::vector<double> x(1, x0 + (xEnd - x0) / 2);
    double h = (xEnd - x0) / 2;

    for (int i = 0; i < 20; i++) {
        int N = static_cast<int>((xEnd - x0) / h) + 1;
        int n = N - 2;
        y = yUpdater(x0, xEnd, alpha, beta, y, x, n, h);
        x = xUpdater(x0, xEnd, alpha, beta, y, x, n, h);
        double error = 1.0;

        std::vector<double> a, b, c, r;
        while (error > accuracy) {
            JCalc(alpha, beta, x, y, n, h, i == 0, a, b, c, r);
            std::vector<double> u(n);
            tridag(a, b, c, r, u);
            y = y + u;
            error = std::sqrt(std::inner_product(u.begin(), u.end(), u.begin(), 0.0));
        }

        std::cout << std::setw(2) << i + 1 << std::setw(15) << y[targetIndex];

        if (i > 0) {
            double diff = prev_y[prev_targetIndex] - y[targetIndex];
            std::cout << std::setw(20) << diff;
            if (i > 1) {
                double diff1 = prev_prev_y[prev_prev_targetIndex] - prev_y[prev_targetIndex];
                double alpha_k = diff1 / diff;
                double richardson = diff / (alpha_k - 1.0);
                std::cout << std::setw(15) << alpha_k
                          << std::setw(15) << richardson
                          << std::setw(15) << std::log2(diff1 / diff);
                if (std::abs(richardson) < accuracy) {
                    std::cout << std::setw(10) << n << std::endl;
                    std::cout << std::endl;
                    return;
                }
            } else {
                std::cout << std::setw(45) << " ";
            }
        } else {
            std::cout << std::setw(65) << " ";
        }
        std::cout << std::setw(10) << n;
        std::cout << std::endl;

        prev_prev_y = prev_y;
        prev_y = y;
        prev_prev_targetIndex = prev_targetIndex;
        prev_targetIndex = targetIndex;
        targetIndex = 2 * targetIndex + 1;
        h /= 2;
    }
    std::cout << std::endl;
}

int main() {
    double x0 = 0.0;
    double xEnd = 2.0;
    double alpha = 0.0;
    double beta = 1.0;
    double accuracy = 1e-6;

    std::cout << "Finite Difference Method" << std::endl;
    finiteDiff(alpha, beta, x0, xEnd, accuracy);
}

