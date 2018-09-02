set terminal pngcairo 
set output 'ex2.png'

set title "ex2 Saving Simple Plots" 
set xlabel "x"
set ylabel "y"

plot [-10:10] sin(x), \
            atan(x) 
