#!/bin/bash

./serial

/usr/bin/mpirun -n 16 ./Broadcast

/usr/bin/mpirun -n 16 ./Scatter

/usr/bin/mpirun -n 16 ./TeamScatter
