// f    : m d^2 x/dt^2 = f(x, v, t)
// x, v : double[steps+1], x[0] = x0, v[0] = v[0]
void rk4(double (*f)(double, double, double), double *x, double *v, double h, int steps){
    for(int i=1;i<=steps;i++){
        double k1[2] = {v[i-1], f(x[i-1], v[i-1], h*(i-1))};
        double k2[2] = {v[i-1] + 0.5*h*k1[1], f(x[i-1] + 0.5*h*k1[0], v[i-1] + 0.5*h*k1[1], h*(i-0.5))};
        double k3[2] = {v[i-1] + 0.5*h*k2[1], f(x[i-1] + 0.5*h*k2[0], v[i-1] + 0.5*h*k2[1], h*(i-0.5))};
        double k4[2] = {v[i-1] + h*k3[1], f(x[i-1] + h*k3[0], v[i-1] + h*k3[1], h*i)};
        x[i] = x[i-1] + h*(k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0])/6.0;
        v[i] = v[i-1] + h*(k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1])/6.0;
    }
}