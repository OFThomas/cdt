
gnuplot << EOF

set terminal gif animate delay 5 loop 0 optimize
set output "rot.gif"

unset surface
set pm3d at s

set label "sin(x)" at screen 0.7, 0.9

n = 100
do for [i=1:n] {
   set view 60, i*360/n
   splot sin(x) notitle
}

set output
EOF
