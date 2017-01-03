set terminal epslatex size 9cm,7cm color colortext
set output "plots/kitaev.tex"

set ylabel '$E/t$'
set xlabel '$\mu/t$'
#set yrange [-0.05: 0.4]
#plot "data.dat" u 1:2 w p, "data.dat" u 1:3 w l, "data.dat" u 1:4 w l
filename(n) = sprintf("dane/data%d.dat", n)
plot for [i=0:49] filename(i) using 1:2 with lines notitle
