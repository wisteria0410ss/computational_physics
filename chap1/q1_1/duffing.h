// constants
    #define PI      3.14159265358979
// params
    #define m       1.0
    #define gamma   0.1
    #define delta   (3.0*PI/2.0)
    #define a       0.25
    #define b       0.5
    #define F0      10.0
    #define omega   (PI/3.0)

// duffing.c
    // ma(t) = F(x, v, t)
    double F(double x, double v, double t);

// rk4.c
    // f    : m d^2 x/dt^2 = f(x, v, t)
    // x, v : double[steps+1], x[0] = x0, v[0] = v[0]
    void rk4(double (*f)(double, double, double), double *x, double *v, double h, int steps);
