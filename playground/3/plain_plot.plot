set terminal epslatex size 9cm,7cm color colortext
set output 'plots/plain/plain.tex'

set linetype 11 lc rgb '#C83B34' lw 2
set linetype 10 lc rgb '#848484'

i = 202
plot './eigenvalues.dat' using 1:i-3 w l ls 10 notitle, './eigenvalues.dat' using 1:i-2 w l ls 10 notitle, './eigenvalues.dat' using 1:i-1 w l ls 11 notitle, './eigenvalues.dat' using 1:i w l ls 11 notitle, './eigenvalues.dat' using 1:i+1 w l ls 10 notitle, './eigenvalues.dat' using 1:i+1 w l ls 10 notitle, './eigenvalues.dat' using 1:i+3 w l ls 10 notitle, './eigenvalues.dat' using 1:i+4 w l ls 10 notitle,
