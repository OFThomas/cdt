gnuplot -e "set key off; set xlabel 'w1'; set ylabel 'w2'; set zlabel 'F(x,y)'; splot 'signalfreq.dat' with linespoints pointtype 7, 'idlerfreq.dat' with linespoints pointtype 6; pause -1"

