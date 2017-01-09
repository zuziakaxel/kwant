set terminal epslatex size 15cm,7cm color colortext
out_name = sprintf("plots/b/b%d.tex", file)
set output out_name

#set ylabel "E/t"
#set xlabel "n"
#set yrange [-0.05: 0.05]
#set xrange [95:105]
#plot "majorana_data.dat" u 1:2 w l
#plot "dane/b/04.dat" u ($1):2 w p

#plot "dane/b/04.dat" u 1:2 w lp, "dane/b/04.dat" u 1:3 w lp, "dane/b/04.dat" u 1:4 w lp, "dane/b/04.dat" u 1:5 w lp
set linetype 10 lc rgb '#0BCC47'
set linetype 11 lc rgb '#C83B34'
filename(n) = sprintf("./dane/b/%d.dat", n)
#splot "data.dat" u 1:2:3 w l
#plot for [i=198:205] filename(file) using 1:i w l ls 10 notitle
i = 201
plot filename(file) using 1:i-3 w l ls 10 notitle, filename(file) using 1:i-2 w l ls 10 notitle, filename(file) using 1:i-1 w l ls 10 notitle, filename(file) using 1:i w l ls 11 notitle, filename(file) using 1:i+1 w l ls 11 notitle, filename(file) using 1:i+1 w l ls 10 notitle, filename(file) using 1:i+3 w l ls 10 notitle, filename(file) using 1:i+4 w l ls 10 notitle,
