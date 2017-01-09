for iter in `seq -30 30`;
do
    gnuplot -e "file=${iter}" main.plot
    # echo $i
done