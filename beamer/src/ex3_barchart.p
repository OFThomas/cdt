set term pngcairo
set output "ex3.png"

unset key
set ylabel "y axis"
set xlabel "x axis"

set boxwidth 1.0
set style fill solid
plot "barchart.dat" with boxes


