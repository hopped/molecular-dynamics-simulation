#set terminal svg size 600,400 dynamic enhanced fname 'arial' fsize 10
set terminal svg size 600,600

set title "Graphite wall"
set view 85,10
set xrange [0:4]
set yrange [0:3]
set zrange [0:3]

set view 82,10
set output "graphite-wall-3D.svg"
splot 'graphite-wall.dat' u 1:2:3 with points pointtype 7 pointsize 1 notitle

set view 0,0
set output "graphite-wall-top.svg"
splot 'graphite-wall.dat' u 1:2:3 with points pointtype 7 pointsize 1 notitle

set view 90,0
set output "graphite-wall-front.svg"
splot 'graphite-wall.dat' u 1:2:3 with points pointtype 7 pointsize 1 notitle

set view 90,90
set output "graphite-wall-side.svg"
splot 'graphite-wall.dat' u 1:2:3 with points pointtype 7 pointsize 1 notitle
