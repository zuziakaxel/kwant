set terminal png
set output "test_plot.png"

set ylabel "E/t"
set xlabel "n"
#set yrange [-5: 5]
#set xrange [95:105]
#plot "majorana_data.dat" u 1:2 w l
plot "eigenvalues.dat" u ($1):2 w p
#plot "eigenvalues.dat" u 1:2 w lp, "eigenvalues.dat" u 1:3 w lp, "eigenvalues.dat" u 1:4 w lp, "eigenvalues.dat" u 1:5 w lp
