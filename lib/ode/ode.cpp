#include <iostream>
#include <functional>

double numerov(std::function<double(double)> f, double t, double x0, double x1, double h){
    double wt = 2.0 * (1.0 + 5.0*h*h*f(t-h)/12.0) * x1 - (1.0 - h*h*f(t-2*h)/12.0) * x0;
    return wt / (1.0 - h*h*f(t)/12.0);
}