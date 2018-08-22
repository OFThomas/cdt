gnuplot << 'EOF'

set terminal gif animate delay 10 loop 0 optimize
set output "g4python.gif"
set label "G^{(4)}" at screen 0.7, 0.9
set hidden3d

unset key

set xlabel 'Squeezing'
set ylabel 'Beamsplitter angle'
set zlabel 'G4'


#set palette model CMY rgbformulae 7,5,15

n = 200
do for [i=1:n] {
   set view 80, i*360/n
    splot 'g4data.txt' using 1:2:3 with linespoints pointtype 7 palette 
}

set output

#pause 30
EOF
