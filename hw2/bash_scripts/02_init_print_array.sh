#!/bin/bash

start=10
end=20

for ((i = 0; i <= $end - $start; ++i)); do my_array[$i]=$(($start + $i)); done

for elem in ${my_array[@]}; do echo $elem; done