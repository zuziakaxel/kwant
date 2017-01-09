for iter in `seq 2 11`;
do
    gnuplot -e "file=${iter}" angle.plot
    # echo $i
done

for iter in `seq 0 10`;
do
    gnuplot -e "file=${iter}" delta.plot
done
