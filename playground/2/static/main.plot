set terminal epslatex size 12cm,7cm color colortext
filename = sprintf("plots/mu_%d.tex", file)
set output filename

set ylabel "E"
set xlabel "k"
set xrange [-3.2: 3.2]
set yrange [-4.1: 4.1]
set samples 10000
#set yrange [-0.05: 0.4]
#plot "data.dat" u 1:2 w p, "data.dat" u 1:3 w l, "data.dat" u 1:4 w l
set linetype 11 lc rgb '#C83B34' lw 2
set linetype 10 lc rgb '#0BCC47'
# filename(n) = sprintf("dane/data%d.dat", n)

t = 1.0
delta = 1.0
mu = file/10.0

f(x) = sqrt( (2*t*cos(x) + mu)**2 + 4*delta**2*sin(x)**2)
g(x) = -sqrt( (2*t*cos(x) + mu)**2 + 4*delta**2*sin(x)**2)
plot f(x) w l ls 11 notitle, g(x) w l ls 11 notitle 
