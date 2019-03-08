#include <functional>

// d^2 x(t) / dt^2 = f(t)x(t) を解く。x0とx1からx2=x(t)を求める。
double numerov(std::function<double(double)> f, double t, double x0, double x1, double h);

// d^2 x(t) / dt^2 = f(t)x(t) を解く。x(t)からx(t+h)を求める。刻み幅hは誤差errに基づき自動適合。
double v_rk4(std::function<double(double)> f, double t, std::pair<double, double> *x0, double err, double *h);