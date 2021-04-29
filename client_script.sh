#!/bin/bash

for k in {0..1}
do
    ./ppDBSCAN_distance  -r 1 -e 400 -d "Lsun_prepared"

    for j in {0..399}
    do
        ./ppDBSCAN_grouping -r 1 -e 400 -i $j -d "Lsun_prepared"
    done
done
