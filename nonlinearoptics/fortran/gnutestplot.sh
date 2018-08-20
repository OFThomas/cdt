FILE=$1
gnuplot << EOF 
plot '$FILE' 
set xlabel 'Freq'
set ylabel 'F()'
set grid
set term dumb
replot 
set term wxt
replot
EOF
