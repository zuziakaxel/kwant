set terminal epslatex size 8cm,7cm color colortext
out_name = sprintf("plots/wave/wave%d.tex", file)
set output out_name

set ylabel '$|u|^2 + |v|^2$'
set xlabel '$n$'
set xrange [0:24]
set yrange [0:0.5]
set linetype 11 lc rgb '#C83B34' lw 1.5
set linetype 12 lc rgb '#0BCC47' lw 1.5
set linetype 14 lc rgb '#353E4E' lw 1.5
set linetype 10 lc rgb '#848484'
# set key outside right center;
filename(n) = sprintf("dane/v/data%d.dat", n)
mu = 4.0*file/200.0
plot filename(file) u 1:2 w l ls 11 notitle
#plot filename(2) u 1:2 w l ls 11 title '$\mu/t = 0.0$',\
#filename(4) u 1:2 w l ls 12 title '$\mu/t = 2.0$',\,\
#filename(8) u 1:2 w l ls 13,
