set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 600, 400 
set output 'ex1.png'
set title "ex1 Simple Plots" 
set title  font ",20" norotate
plot [-10:10] sin(x),atan(x),cos(atan(x))
