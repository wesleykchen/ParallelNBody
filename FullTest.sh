#!/bin/bash

mpirun -n 16 ./symmetric -c 4 data/oneD_source_1000 data/oneD_charge_1000

mpirun -n 16 ./teamscatter -c 4 data/oneD_source_1000 data/oneD_charge_1000

./serial data/oneD_source_1000 data/oneD_charge_1000

mpirun -n 16 ./symmetric -c 4 data/oneD_source_10000 data/oneD_charge_10000

mpirun -n 16 ./teamscatter -c 4 data/oneD_source_10000 data/oneD_charge_10000

./serial data/oneD_source_10000 data/oneD_charge_10000

mpirun -n 16 ./symmetric -c 4 data/oneD_source_100000 data/oneD_charge_100000

mpirun -n 16 ./teamscatter -c 4 data/oneD_source_100000 data/oneD_charge_100000

./serial data/oneD_source_100000 data/oneD_charge_100000