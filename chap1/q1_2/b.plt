set terminal png
set output "b.png"

a = 1.0
b = 1.0
f(x) = a*x+b
fit f(x) "b.dat" using (log($1)):(log($2)) via a,b

set xlabel "R_g"
set ylabel "N"
set logscale x
set logscale y
plot "b.dat" title "data", exp(b)*x**a title "fitting"