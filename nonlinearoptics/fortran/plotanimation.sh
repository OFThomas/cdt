gnuplot << 'EOF'

set terminal gif animate delay 5 loop 0 optimize
set output "jsa.gif"
set label "JSA" at screen 0.7, 0.9

unset key
set palette model CMY rgbformulae 7,5,15
n = 100
do for [i=1:n] {
   set view 60, i*360/n
    splot 'fplotw1w2.dat' using 1:2:($3 <=0.05 ? NaN : $3) with linespoints palette, 'fplotw3w4.dat' using 1:2:($3 <=0.05 ? NaN : $3) with linespoints palette
}

set output

#pause 30
EOF
