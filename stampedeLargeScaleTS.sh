#!/bin/bash
#SBATCH --ntasks 4096         #Number of processes
#SBATCH -t 01:00:00              #Runtime in minutes
#SBATCH -p normal   	      #Partition to submit to
#SBATCH -o data/stampedeTS256k.out     	      #File to which standard out will be written
#SBATCH -e data/stampedeTS256k.err      	      #File to which standard err will be written

module swap intel gcc/4.7.1

make clean

make teamscatter XFLAGS='-DP2P_NUM_THREADS=0 -DP2P_DECAY_ITERATOR=0'

#
# Execute the run
#

mpirun -np 4096  ./teamscatter 256000  -c 1 #-nocheck

mpirun -np 4096  ./teamscatter 256000  -c 2 #-nocheck

mpirun -np 4096  ./teamscatter 256000  -c 4 #-nocheck

mpirun -np 4096  ./teamscatter 256000  -c 8 #-nocheck

mpirun -np 4096  ./teamscatter 256000  -c 16 #-nocheck

mpirun -np 4096  ./teamscatter 256000  -c 32 #-nocheck

mpirun -np 4096  ./teamscatter 256000  -c 64 #-nocheck