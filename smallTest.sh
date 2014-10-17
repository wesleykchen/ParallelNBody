#!/bin/bash
#SBATCH --ntasks 64         #Number of processes
#SBATCH -t 300                #Runtime in minutes
#SBATCH -p normal   	      #Partition to submit to

#SBATCH --mem-per-cpu=200     #Memory per cpu in MB (see also --mem)
#SBATCH -o data/largeScale256k.out     	      #File to which standard out will be written
#SBATCH -e data/largeScale256k.err      	      #File to which standard err will be written

module swap intel gcc/4.7.1

make clean

make symmetric XFLAGS='-DP2P_NUM_THREADS=0 -DP2P_DECAY_ITERATOR=0'

#
# Execute the run
#

mpirun -np 64  ./symmetric 25600  -c 4 #-nocheck