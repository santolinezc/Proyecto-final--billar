reset
unset key
set size ratio -1
set parametric
set trange [0:1]
set term gif animate delay 2
unset ytics
unset xtics
unset border
set output 'Billar_trajectories.gif'
counter = 0
M = 1
alpha = 0.01
R = 3
do for[ii = 1:999997:3 ]{if( ii < 999997/50 ){plot for [n = 0:M] 'Data.txt' using 2+2*n:3+2*n every ::1::ii w l ls n lt 8 dashtype 2, R*cos(pi*t), alpha + R*sin(pi*t) lt 2 lw 2, R*cos(pi*t), -(alpha + R*sin(pi*t)) lt 2 lw 2, R, -alpha + 2*alpha*t lt 2 lw 2, -R, -alpha + 2*alpha*t lt 2 lw 2, for [n = 0:M] 'Data.txt'  using 2+2*n:3+2*n every ::ii::ii w p ls 6 ps 2 pt 7} else {plot for [n = 0:M] 'Data.txt' using 2+2*n:3+2*n every ::(1+counter)::ii w l ls n lt 8, for [n = 0:M] 'Data.txt'  using 2+2*n:3+2*n every::ii::ii w p ls 6 ps 2 pt 7, R*cos(pi*t), alpha + R*sin(pi*t) lt 2 lw 2, R*cos(pi*t), -(alpha + R*sin(pi*t)) lt 2 lw 2, R, -alpha + 2*alpha*t lt 2 lw 2, -R, -alpha + 2*alpha*t lt 2 lw 2; counter = counter + 1}}
