set terminal png
set output "d.png"

set xlabel "b"
set ylabel "N"

a=2.0
b=1.0
f(x) = a*x+b
fit f(x) "d.dat" using (log($1)):(log($2)) via a,b

set logscale x
set logscale y
set xr [6.0/128.0/1.1:6.0/2.0*1.1]

plot "d.dat" title "data", exp(b)*x**a title "fitting"