set terminal png
set output "b.png"

set xlabel "x"
set ylabel "p"

plot "b.dat" using 2:3 pt 0 notitle #title "Strange Attractor"