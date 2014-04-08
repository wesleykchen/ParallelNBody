#!/bin/bash

mpirun -n 8 ./symmetric 1000 -c 2

mpirun -n 8 ./teamscatter 1000 -c 2

./serial 1000

mpirun -n 8 ./symmetric 10000 -c 2

mpirun -n 8 ./teamscatter 10000 -c 2

mpirun -n 16 ./symmetric 10000 -c 2

mpirun -n 16 ./teamscatter 10000 -c 2

mpirun -n 16 ./symmetric 10000 -c 4

mpirun -n 16 ./teamscatter 10000 -c 4

./serial 10000 -nocheck

mpirun -n 8 ./symmetric 100000 -c 2 -nocheck

mpirun -n 8 ./teamscatter 100000 -c 2 -nocheck

mpirun -n 16 ./symmetric 100000 -c 2 -nocheck

mpirun -n 16 ./teamscatter 100000 -c 2 -nocheck

mpirun -n 16 ./symmetric 100000 -c 4 -nocheck

mpirun -n 16 ./teamscatter 100000 -c 4 -nocheck

./serial 100000 -nocheck


