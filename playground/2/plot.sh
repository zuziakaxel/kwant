#!/bin/bash
for iter in `seq 1 200`;
do
    gnuplot -e "file=${iter}" wave.plot
    # echo $i
done
