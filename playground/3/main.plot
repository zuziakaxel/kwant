set terminal png
set output "test_plot.png"

#set ylabel "E/t"
#set xlabel "n"
#set yrange [-0.05: 0.05]
#set xrange [95:105]
#plot "majorana_data.dat" u 1:2 w l
#plot "eigenvalues.dat" u ($1):2 w p

#plot "eigenvalues.dat" u 1:2 w lp, "eigenvalues.dat" u 1:3 w lp, "eigenvalues.dat" u 1:4 w lp, "eigenvalues.dat" u 1:5 w lp
set linetype 10 lc rgb '#0BCC47'
set linetype 11 lc rgb '#C83B34'
#splot "data.dat" u 1:2:3 w l
#plot for [i=198:205] './eigenvalues.dat' using 1:i w l ls 10 notitle
i = 201
plot './eigenvalues.dat' using 1:i-3 w l ls 10 notitle, './eigenvalues.dat' using 1:i-2 w l ls 10 notitle, './eigenvalues.dat' using 1:i-1 w l ls 10 notitle, './eigenvalues.dat' using 1:i w l ls 11 notitle, './eigenvalues.dat' using 1:i+1 w l ls 11 notitle, './eigenvalues.dat' using 1:i+1 w l ls 10 notitle, './eigenvalues.dat' using 1:i+3 w l ls 10 notitle, './eigenvalues.dat' using 1:i+4 w l ls 10 notitle,
