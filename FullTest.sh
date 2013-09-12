#!/bin/bash

./serial

mpirun -n 16 ./broadcast

mpirun -n 16 ./scatter

mpirun -n 16 ./teamscatter
