#!/bin/bash

chmod +x mpi_cellular_automata.py

for i in $(seq 1 $1)
do
    mpirun -n $i ./mpi_cellular_automata.py --n-epochs 12000 --size 720 --rule-id 110 --periodic
done

echo "Done!"