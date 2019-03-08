set terminal png
set output "differential_cross_section.png"

set xrange [0:pi]
set xlabel "theta [rad]"
set ylabel "Differential cross section"

plot "differential_cross_section.dat" with lines title "Differential cross section"