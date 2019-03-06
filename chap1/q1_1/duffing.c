#include <math.h>
#include "duffing.h"

// ma(t) = F(x, v, t)
double F(double x, double v, double t){
    return -gamma*v + 2*a*x - 4*b*pow(x, 3.0) + F0*cos(omega*t + delta);
}