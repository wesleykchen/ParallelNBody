#!/bin/bash
#SBATCH --ntasks 2048         #Number of processes
#SBATCH -t 01:00:00                #Runtime in minutes
#SBATCH -p normal   	      #Partition to submit to
#SBATCH -o data/a_stampedefull256k_wait.out     	      #File to which standard out will be written
#SBATCH -e data/a_stampedefull256k_wait.err      	      #File to which standard err will be written

module load intel gcc/4.7.1

make clean

make symmetric XFLAGS='-DP2P_NUM_THREADS=0 -DP2P_DECAY_ITERATOR=0'

#
# Execute the run
#
ibrun -n 1  -o 0 ./symmetric 256000  -c 1 #-nocheck

wait

ibrun -n 2  -o 0 ./symmetric 256000  -c 1 #-nocheck

wait

ibrun -n 4  -o 0 ./symmetric 256000  -c 1 #-nocheck

wait

ibrun -n 4  -o 0 ./symmetric 256000  -c 2 #-nocheck

wait

ibrun -n 8  -o 0 ./symmetric 256000  -c 1 #-nocheck

wait

ibrun -n 8  -o 0 ./symmetric 256000  -c 2 #-nocheck

wait

ibrun -n 16  -o 0 ./symmetric 256000  -c 1 #-nocheck

wait

ibrun -n 16  -o 0 ./symmetric 256000  -c 2 #-nocheck

wait

ibrun -n 16  -o 0 ./symmetric 256000  -c 4 #-nocheck

wait

ibrun -n 32  -o 0 ./symmetric 256000  -c 1 #-nocheck

wait

ibrun -n 32  -o 0 ./symmetric 256000  -c 2 #-nocheck

wait

ibrun -n 32  -o 0 ./symmetric 256000  -c 4 #-nocheck

wait

ibrun -n 64  -o 0 ./symmetric 256000  -c 1 #-nocheck

wait

ibrun -n 64  -o 0 ./symmetric 256000  -c 2 #-nocheck

wait

ibrun -n 64  -o 0 ./symmetric 256000  -c 4 #-nocheck

wait

ibrun -n 64  -o 0 ./symmetric 256000  -c 8 #-nocheck

wait

ibrun -n 128  -o 0 ./symmetric 256000  -c 1 #-nocheck

wait

ibrun -n 128  -o 0 ./symmetric 256000  -c 2 #-nocheck

wait

ibrun -n 128  -o 0 ./symmetric 256000  -c 4 #-nocheck

wait

ibrun -n 128  -o 0 ./symmetric 256000  -c 8 #-nocheck

wait

ibrun -n 256  -o 0 ./symmetric 256000  -c 1 #-nocheck

wait

ibrun -n 256  -o 0 ./symmetric 256000  -c 2 #-nocheck

wait

ibrun -n 256  -o 0 ./symmetric 256000  -c 4 #-nocheck

wait

ibrun -n 256  -o 0 ./symmetric 256000  -c 8 #-nocheck

wait

ibrun -n 256  -o 0 ./symmetric 256000  -c 16 #-nocheck

wait

ibrun -n 512  -o 0 ./symmetric 256000  -c 1 #-nocheck

wait

ibrun -n 512  -o 0 ./symmetric 256000  -c 2 #-nocheck

wait

ibrun -n 512  -o 0 ./symmetric 256000  -c 4 #-nocheck

wait

ibrun -n 512  -o 0 ./symmetric 256000  -c 8 #-nocheck

wait

ibrun -n 512  -o 0 ./symmetric 256000  -c 16 #-nocheck

wait

ibrun -n 1024  -o 0 ./symmetric 256000  -c 1 #-nocheck

wait

ibrun -n 1024  -o 0 ./symmetric 256000  -c 2 #-nocheck

wait

ibrun -n 1024  -o 0 ./symmetric 256000  -c 4 #-nocheck

wait

ibrun -n 1024  -o 0 ./symmetric 256000  -c 8 #-nocheck

wait

ibrun -n 1024  -o 0 ./symmetric 256000  -c 16 #-nocheck

wait

ibrun -n 1024  -o 0 ./symmetric 256000  -c 32 #-nocheck

wait

ibrun -n 2048  -o 0 ./symmetric 256000  -c 1 #-nocheck

wait

ibrun -n 2048  -o 0 ./symmetric 256000  -c 2 #-nocheck

wait

ibrun -n 2048  -o 0 ./symmetric 256000  -c 4 #-nocheck

wait

ibrun -n 2048  -o 0 ./symmetric 256000  -c 8 #-nocheck

wait

ibrun -n 2048  -o 0 ./symmetric 256000  -c 16 #-nocheck

wait

ibrun -n 2048  -o 0 ./symmetric 256000  -c 32 #-nocheck

wait
