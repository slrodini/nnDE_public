#!/usr/bin/gnuplot
set terminal postscript eps enhanced color

set lmargin 13
set rmargin 3
set tmargin 3
set bmargin 6

#set for [i=1:2] linetype i dt i

set style line 1 dt 1 linecolor rgbcolor "#4B0082" linewidt 5
set style line 3 dt 3 linecolor rgbcolor "#0000FF" linewidt 5
set style line 2 dt 5 linecolor rgbcolor "#FF0000" linewidt 5
set style line 4 dt 1 linecolor rgbcolor "#FFFFFF" linewidt 1
set style line 5 dt 1 linecolor rgbcolor "#000000" linewidt 5
set style line 6 dt 2 linecolor rgbcolor "#000000" linewidt 5

set ylabel "log(t/t_0)" offset -2,-0.5  font ",28" rotate by 90 tc rgbcolor "#000000"
set xlabel "N_0" font ",28" offset 0,-1
set ytics font ",24"
set format y "%1.1f"
#set ytics add ("0" 0)
set xtics font ",24" offset 0,-0.5

set xrange[:60]
set yrange[-12:-1]
set key top left notitle spacing 1.5 font ",24"

f(x) = a*x+b*sqrt(x)+c
g(x) = d*x+e*sqrt(x)+p
h(x) = q*x+r*sqrt(x)+s

set output "time_nI.eps"
set xlabel "N_0" font ",28" offset 5,-1
fit [x=1:60]  f(x) 'time_nI.dat' u 1:(log($2)) via a,b,c
fit [x=1:60]  g(x) 'time_nI.dat' u 1:(log($3)) via d,e,p
fit [x=1:60]  h(x) 'time_nI.dat' u 1:(log($4)) via q,r,s
plot 'time_nI.dat' u 1:(log($2)) pointsize 2 pointtype 12 linecolor rgbcolor "#4B8200" title 'full Hessian',\
        f(x) w l ls 5 notitle,\
     'time_nI.dat' u 1:(log($3)) pointsize 2 pointtype 10 linecolor rgbcolor "#4B0082" title 'diagonal Hessian',\
        g(x) w l ls 6 notitle,\
     'time_nI.dat' u 1:(log($4)) pointsize 2 pointtype 8 linecolor rgbcolor "#004B82" title 'first derivatives',\
        h(x) w l ls 6 notitle 


unset yrange
set output "time_nO.eps"
set key bottom right notitle spacing 1.5 font ",24"


set xlabel "N_L" font ",28" offset 0,-1

fit [x=1:60]  f(x) 'time_nO.dat' u 1:(log($2)) via a,b,c
fit [x=1:60]  g(x) 'time_nO.dat' u 1:(log($3)) via d,e,p
fit [x=1:60]  h(x) 'time_nO.dat' u 1:(log($4)) via q,r,s
plot 'time_nO.dat' u 1:(log($2)) pointsize 2 pointtype 12 linecolor rgbcolor "#4B8200" title 'full Hessian',\
        f(x) w l ls 5 notitle,\
     'time_nO.dat' u 1:(log($3)) pointsize 2 pointtype 10 linecolor rgbcolor "#4B0082" title 'diagonal Hessian',\
        g(x) w l ls 6 notitle,\
     'time_nO.dat' u 1:(log($4)) pointsize 2 pointtype 8 linecolor rgbcolor "#004B82" title 'first derivatives',\
        h(x) w l ls 6 notitle

