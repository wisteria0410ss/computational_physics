#include <iostream>
#include <utility>
#include <cmath>
#include <functional>

double numerov(std::function<double(double)> f, double t, double x0, double x1, double h){
    double wt = 2.0 * (1.0 + 5.0*h*h*f(t-h)/12.0) * x1 - (1.0 - h*h*f(t-2*h)/12.0) * x0;
    return wt / (1.0 - h*h*f(t)/12.0);
}

std::pair<double, double> operator*(double a, std::pair<double, double> x){
    return {a*x.first, a*x.second};
}
std::pair<double, double> operator*(std::pair<double, double> x, double a){
    return {a*x.first, a*x.second};
}
std::pair<double, double> operator/(std::pair<double, double> x, double a){
    return {x.first/a, x.second/a};
}
std::pair<double, double> operator+(std::pair<double, double> x, std::pair<double, double> y){
    return {x.first+y.first, x.second+y.second};
}

std::pair<double, double> rk4(std::function<double(double)> f, double t, std::pair<double, double> x0, double h){
    std::pair<double, double> k1 = {h * f(t) * x0.second, h * x0.first};
    std::pair<double, double> k2 = {h * f(t + h/2.0) * (x0.second + k1.first/2.0), h * (x0.first + k1.second/2.0)};
    std::pair<double, double> k3 = {h * f(t + h/2.0) * (x0.second + k2.first/2.0), h * (x0.first + k2.second/2.0)};
    std::pair<double, double> k4 = {h * f(t + h) * (x0.second + k3.first), h * (x0.first + k3.second)};
    return x0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
}

double v_rk4(std::function<double(double)> f, double t, std::pair<double, double> *prev, double err, double *h){
    auto x0 = *prev;
    double next_1 =rk4(f, t, x0, *h).second;
    auto x1 = rk4(f, t, x0, *h/2.0);
    double next_2 = rk4(f, t + *h/2.0, x1, *h/2.0).second;

    double Delta = std::abs(next_1 - next_2);
    *h = 15.0/16.0 * (*h) * std::pow(err/Delta, 0.2);
    *prev = rk4(f, t, x0, *h);
    return prev->second;
}