// d^2 x(t) / dt^2 = f(t)x(t) を解く。x0とx1からx2=x(t)を求める。
double numerov(double (*f)(double), double t, double n0, double n1, double h);