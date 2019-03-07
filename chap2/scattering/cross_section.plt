set terminal png
set output "cross_section.png"

set xlabel "Energy [meV]"
set ylabel "Cross section [\\AA^2]"

plot "cross_section.dat" with lines title "Cross section"