set term pngcairo
set output "ex4.png"

set multiplot layout 1,2 

set yrange [-0.3:0.35]
set xrange [0:100]

set key box opaque

set ylabel 'y axis'
set xlabel 'x axis'

set style line 1 lw 3 lc 1
set style line 2 lw 3 lc 2
set style line 3 lw 3 lc 6

set style line 4 lw 3 lc 3
set style line 5 lw 3 lc 5
set style line 6 lw 3 lc 2
unset key

#plot for [IDX=0:*] 'multiplot1.dat' i IDX using 2:3 with lines ls (IDX+1)

#plot for [IDX=0:*] 'multiplot2.dat' i IDX using 2:3 with lines ls (IDX+4)

plot 'multiplot1.dat' using 2:3 with lines ls 1
plot 'multiplot2.dat' using 2:3 with lines ls 2

unset multiplot

