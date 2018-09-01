set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 600, 400 
set output 'ex2.png'

set title "ex2 Saving Simple Plots" 
set xlabel "x"
set ylabel "y"

plot [-10:10] sin(x), \
            atan(x) 
