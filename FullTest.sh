#!/bin/bash

mpirun -n 16 ./symmetric -c 2 data/oneD_source_10000 data/oneD_charge_10000

mpirun -n 16 ./symmetric -c 2 data/oneD_source_100000 data/oneD_charge_100000

mpirun -n 8 ./symmetric -c 2 data/oneD_source_100000 data/oneD_charge_10000

mpirun -n 8 ./symmetric -c 2 data/oneD_source_100000 data/oneD_charge_100000

