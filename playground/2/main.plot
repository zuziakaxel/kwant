set terminal png
set output "test_plot.png"

set ylabel "E/t"
set xlabel "mu/t"
plot "data.dat" u 1:2 w p
