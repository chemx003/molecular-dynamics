#set terminal gif animate delay 100 
#set output 'nematic.gif'
set terminal qt size 900, 900;
set view 85, 0
set xrange [-5:20]
set yrange [-5:20]
set zrange [-5:20]
do for [i=2:999] {print i; splot sprintf('vector%i.data', 1*i) using 1:2:3:4:5:6 with vectors nohead; pause 0.01}
pause -1
