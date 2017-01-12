#!/bin/bash

for i in `seq 1 3`
do
    sed "s/indk\ =\ 1/indk\ =\ $i/g" 3D_MSE_bc_Ksolver.edp > tmp.edp
    FreeFem++-nw tmp.edp >> /dev/null
    rm tmp.edp
    echo "Done K problem indf=$i"
    sed "s/indk\ =\ 1/indk\ =\ $i/g" 3D_MSE_bc_Lsolver.edp > tmp.edp
    FreeFem++-nw tmp.edp >> /dev/null
    echo "Done L problem indk=$i, indl=3"
    rm tmp.edp
done
