set terminal png
set output "test_plot.png"

set ylabel "E/t"
set xlabel "mu/t"
set yrange [-0.05: 0.4]
plot "majorana_data.dat" u 1:2 w l
