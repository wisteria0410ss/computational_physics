set terminal png
set output "a.png"

set xlabel "t"
set ylabel "x"

plot "a.dat" using 1:2 with lines title "Duffing"