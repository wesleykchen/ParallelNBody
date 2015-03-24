#!/bin/bash
#SBATCH --ntasks 256         #Number of processes
#SBATCH -t 00:15:00                #Runtime in minutes
#SBATCH -p normal   	      #Partition to submit to
#SBATCH -o data/a_stampede256kp256.out     	      #File to which standard out will be written
#SBATCH -e data/a_stampede256kp256.err      	      #File to which standard err will be written

module load intel gcc/4.7.1

make clean

make symmetric XFLAGS='-DP2P_NUM_THREADS=0 -DP2P_DECAY_ITERATOR=0'

#
# Execute the run
#

ibrun -n 256 -o 0 ./symmetric 256000  -c 1 #-nocheck

ibrun -n 256 -o 0 ./symmetric 256000  -c 2 #-nocheck

ibrun -n 256 -o 0 ./symmetric 256000  -c 4 #-nocheck

ibrun -n 256 -o 0 ./symmetric 256000  -c 8 #-nocheck

ibrun -n 256 -o 0 ./symmetric 256000  -c 16 #-nocheck

ibrun -n 256 -o 0 ./symmetric 256000  -c 32 #-nocheck

