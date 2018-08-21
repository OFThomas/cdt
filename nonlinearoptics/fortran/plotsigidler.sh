#gnuplot -e "set xlabel 'w1'; set ylabel 'w2'; set zlabel 'F(x,y)'; splot 'signalfreq1.dat' with linespoints pointtype 7, 'signalfreq2.dat' with linespoints pointtype 6; pause -1"

gnuplot << EOF
set term png
set output "single_sig_idler12.png"

set multiplot layout 1,2 

#set format y "%.1f"
set yrange [0:0.35]

<<<<<<< HEAD
=======
#set format y "%.1f"
>>>>>>> 732d5106599e5bbe2a52e471720567f5638835e3
set key box opaque

set ylabel 'Amplitude'
set xlabel 'Frequency'

set style line 1 lc 1
set style line 2 lc 3

unset key

plot 'signalfreq1.dat' using 2:3 with lines ls 1
plot 'idlerfreq1.dat' using 2:3 with lines ls 2

plot 'signalfreq2.dat' using 2:3 with lines ls 1
plot 'idlerfreq2.dat' using 2:3 with lines ls 2

unset ylabel
unset ytics

unset multiplot
#pause 1
EOF

<<<<<<< HEAD
# show plots while waiting for animation
display single_sig_idler12.png & 

gnuplot << 'EOF'

set terminal gif animate delay 10 loop 0 optimize
set output "jsa12.gif"
set label "JSA" at screen 0.7, 0.9

unset key
set hidden3d
#set palette model CMY rgbformulae 7,5,15

set xlabel "w_{s}"
set ylabel "w_{i}"
set zlabel "Amplitude"

set ticslevel 0

n = 100
do for [i=1:n] {
   set view 60, i*360/n
    splot 'fplotw1w2.dat' using 1:2:($3 <=0.000000000000000000000000000000000000000000000000005 ? NaN : $3) with linespoints palette pointtype 7, 'fplotw3w4.dat' using 1:2:3 with linespoints palette pointtype 5
}

set output

#pause 30
EOF

animate jsa12.gif 
=======

>>>>>>> 732d5106599e5bbe2a52e471720567f5638835e3
