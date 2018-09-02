set term pngcairo
set output "ex3.png"

unset key
set ylabel "Occupation"
set xlabel "Schimdt mode number"

set boxwidth 1.0
set style fill solid
plot "schmidtout.dat" with boxes


