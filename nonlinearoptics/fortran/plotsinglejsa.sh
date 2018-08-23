
gnuplot << EOF
set term png
set output "single_sig_idler1.png"

set multiplot layout 1,2 

#set format y "%.1f"
set yrange [-0.3:0.35]
#set xrange [20:105]

set key box opaque

set ylabel 'Amplitude'
set xlabel 'Frequency'

set style line 1 lc 1
set style line 2 lc 3

unset key

plot for [IDX=0:*] 'signalfreq1.dat' i IDX using 2:3 with lines linewidth 3

plot for [IDX=0:*] 'idlerfreq1.dat' i IDX using 2:3 with lines linewidth 3

#plot 'signalfreq1.dat' using 2:3 with lines ls 1
#plot 'idlerfreq1.dat' using 2:3 with lines ls 2

unset ylabel
unset ytics

unset multiplot
#pause 1
EOF

# show plots while waiting for animation
display single_sig_idler1.png & 

#plot schmidt modes

gnuplot << EOF
set term png
set output "schmidtmodesocc.png"

unset key
set ylabel "Occupation"
set xlabel "Schimdt mode number"

set boxwidth 1.0
set style fill solid
plot "schmidtout.dat" with boxes

EOF

# ahow occupation
display schmidtmodesocc.png &

gnuplot << 'EOF'

set terminal gif animate delay 10 loop 0 optimize
set output "jsa1.gif"
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
    splot 'fplotw1w2.dat' using 1:2:($3 <=0.000000000000000000000000000000000000000000000000005 ? NaN : $3) with linespoints palette pointtype 7
}

set output

#pause 30
EOF

animate jsa1.gif 
