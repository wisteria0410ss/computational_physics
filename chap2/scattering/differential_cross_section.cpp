#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>
#include <gsl/gsl_specfunc.h>
#include "../../headers/ode.h"

// lに関する和を打ち切る閾値
const double eps = 1e-9;
// 積分の刻み幅
const double dr = 5e-4;

const double sigma = 3.57;      // Å
const double epsilon = 5.9;     // meV
const double alpha = 6.12/pow(sigma, 2.0);
const double C = sqrt(epsilon*alpha/25.0)*pow(sigma, 6.0);

double V(double r){
    return epsilon * (pow(sigma/r, 12.0) - 2.0*pow(sigma/r, 6.0));
}

double F(int l, double r, double E){
    if(l==0) return (V(r) - E) * alpha;
    else return (V(r) + l*(l+1)/(alpha*r*r) - E) * alpha;
}

void get_wf(double h, int l, double E, double rmax, double *r1, double *r2, double *u_r1, double *u_r2){
    int step = 0;
    double lambda = 2*M_PI/sqrt(alpha*E);

    const double r_min = 0.6 * sigma;
    double prev[2] = {exp(-C*pow(r_min, -5)), exp(-C*pow(r_min + h, -5))};
    for(double r = r_min+h*2; ; r+=h){
        double u = numerov(std::bind(F, l, std::placeholders::_1, E), r, prev[0], prev[1], h);
        //std::cerr << r << "\t" << u << std::endl;
        if(step == 0 && r > rmax){
            step  = 1;
            *r1   = r;
            *u_r1 = u;
        }else if(step == 1 && r >= *r1 + lambda*0.4){
            *r2   = r;
            *u_r2 = u;
            break;
        }
        prev[0] = prev[1];
        prev[1] = u;
    }

    return;
}

double delta(int l, double E, double r1, double r2, double u_r1, double u_r2){
    const double k = sqrt(alpha * E);
    double K = (r1*u_r2) / (r2*u_r1);
    return atan((K*gsl_sf_bessel_jl(l, k*r1) - gsl_sf_bessel_jl(l, k*r2)) / (K*gsl_sf_bessel_yl(l, k*r1) - gsl_sf_bessel_yl(l, k*r2)));
}

double differential_cross_section(double E, double theta){
    double dcs_re = 0.0, dcs_im = 0.0;
    for(int l=0;;l++){
        double r1, r2, u_r1, u_r2;
        get_wf(dr, l, E, 5.0*sigma, &r1, &r2, &u_r1, &u_r2);
        double dl = delta(l, E, r1, r2, u_r1, u_r2);
        double dcs_l = (2*l+1) * sin(dl) * gsl_sf_legendre_Pl(l, cos(theta));
        dcs_re += cos(dl) * dcs_l;
        dcs_im += sin(dl) * dcs_l;
        if(abs(dcs_l) < eps) break;
    }
    return (dcs_re*dcs_re + dcs_im*dcs_im)/sqrt(alpha * E);
}

int main(){
    const double th_min = 0.0, th_max = M_PI;
    const double dth = 0.01;
    double E = 1.0;
    const int N = (th_max - th_min)/dth + 1;
    double theta[N], dcs[N];
    int cnt = 0;
    #pragma omp parallel for
    for(int i=0;i<N;i++){
        double th = th_min + dth*i;
        theta[i] = th;
        dcs[i] = differential_cross_section(E, th);
        #pragma omp critical
        {
            cnt++;
            std::cerr << "\r" << cnt << "/" << N;
        }

    }
    std::cerr << std::endl;

    for(int i=0;i<N;i++){
        std::cout << theta[i] << "\t" << dcs[i] << std::endl;
    }
/*
    int l = 3;
    double E = 1.0;
    const double k = sqrt(alpha * E);

    double r1, r2, u_r1, u_r2;
    get_wf(dr, l, E, 10.0*sigma, &r1, &r2, &u_r1, &u_r2);
    double dl = delta(l, E, r1, r2, u_r1, u_r2);
    for(double r=sigma;r<=r2;r+=dr){
        std::cout << r << "\t" << k*r*(cos(dl)*gsl_sf_bessel_jl(l, k*r) - sin(dl)*gsl_sf_bessel_yl(l, k*r)) << std::endl;
        //std::cout << r << "\t" << sin(k*r - l*M_PI/2.0 + dl) << std::endl;
    }    
*/
    return 0;
}