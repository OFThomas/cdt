set terminal pngcairo 
set output "ex5.png"

unset key
set hidden3d
set palette model CMY rgbformulae 7,5,15

set xlabel "x axis"
set ylabel "y axis"
set zlabel "z axis"

splot 'splot.dat' using 1:2:3 with linespoints palette


