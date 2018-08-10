# basic plot
#gnuplot -e " splot 'g4data.txt' ; pause -1"

# very colourful but lags
gnuplot -e "set key off;set xlabel 'Squeezing'; set ylabel 'Beamsplitter angle'; set zlabel 'G4'; splot 'g4data.txt' with linespoints palette pointtype 7 ; pause -1"

# 2d top view
#gnuplot -e "set key off; set view map; set size ratio 0.9; set object 1 rect from graph 0, graph 0 to graph 1, graph 1 back; set object 1 rect fc rgb 'black' fillstyle solid 1.0; splot 'g4data.txt' with points pointtype 5 pointsize 1 palette linewidth 30; pause -1 "
