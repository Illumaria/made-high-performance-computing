#!/bin/bash

FLAG_ARR=("-g" "-O3")
N_ARR=(500 512 1000 1024 2000 2048)

for flag in ${FLAG_ARR[@]}; do
    echo "Flag: $flag"
    make -s OPTIMIZE=$flag
    for N in ${N_ARR[@]}; do
        make -s run N=$N
    done
done

echo Done!