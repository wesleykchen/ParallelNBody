#!/bin/bash
#SBATCH --ntasks 2048         #Number of processes
#SBATCH -t 01:00:00              #Runtime in minutes
#SBATCH -p normal   	      #Partition to submit to
#SBATCH -o data/a_stampedeTS256k_t2.out     	      #File to which standard out will be written
#SBATCH -e data/a_stampedeTS256k_t2.err      	      #File to which standard err will be written

module swap intel gcc/4.7.1
	
make clean

make teamscatter XFLAGS='-DP2P_NUM_THREADS=0 -DP2P_DECAY_ITERATOR=0'

#
# Execute the run
#

ibrun -n 2048  -o 0 ./teamscatter 256000  -c 1 #-nocheck

ibrun -n 2048  -o 0 ./teamscatter 256000  -c 2 #-nocheck

ibrun -n 2048  -o 0 ./teamscatter 256000  -c 4 #-nocheck

ibrun -n 2048  -o 0 ./teamscatter 256000  -c 8 #-nocheck

ibrun -n 2048  -o 0 ./teamscatter 256000  -c 16 #-nocheck

ibrun -n 2048  -o 0 ./teamscatter 256000  -c 32 #-nocheck