#include <iostream>
#include <cmath>
#include "../../headers/ode.h"

const double sigma = 3.57;      // Å
const double epsilon = 5.9;     // meV
const double alpha = 6.12/std::pow(sigma, 2.0);
const double C = std::sqrt(epsilon*alpha/25.0)*std::pow(sigma, 6.0);

double V(double r){
    return epsilon * (std::pow(sigma/r, 12.0) - 2.0*std::pow(sigma/r, 6.0));
}

double F(int l, double r, double E){
    if(l==0) return (V(r) - E) * alpha;
    else return (V(r) + l*(l+1)/(alpha*r*r) - E) * alpha;
}

void wf_numerov(double h, int l, double E, double rmax){
    double lambda = 2*M_PI/sqrt(alpha*E);

    const double r_min = 0.6 * sigma;
    double prev[2] = {exp(-C*pow(r_min - h, -5)), exp(-C*pow(r_min, -5))};
    for(double r = r_min+h; ; r+=h){
        double u = numerov(std::bind(F, l, std::placeholders::_1, E), r, prev[0], prev[1], h);
        std::cerr << r << "\t" << u << std::endl;
        if(r > rmax){
            break;
        }
        prev[0] = prev[1];
        prev[1] = u;
    }

    return;
}

void wf_rk4(double h, int l, double E, double rmax){
    double lambda = 2*M_PI/sqrt(alpha*E);

    const double r_min = 0.6 * sigma;
    std::pair<double, double> prev = {5*C*std::pow(r_min, -6.0)*std::exp(-C*std::pow(r_min, -5)), std::exp(-C*std::pow(r_min, -5))};
    double r = r_min;
    while(r <= rmax){
        double h = 1e-3;
        double u = v_rk4(std::bind(F, l, std::placeholders::_1, E), r, &prev, 1e-12, &h);
        r += h;
        std::cout << r << "\t" << u << std::endl;
    }
    std::cerr << std::endl;

    return;
}

int main(){
    wf_numerov(1e-3, 0, 1.0, 5*sigma);
    wf_rk4(1e-6, 0, 1.0, 5*sigma);
}