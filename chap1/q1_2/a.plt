set terminal png size 640,640
set output "a.png"

unset colorbox
unset xtics
unset ytics

unset border
set margins 1,1,1,1
set size square
plot "a.dat" matrix with image