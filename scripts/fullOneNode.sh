#!/bin/bash
#SBATCH --ntasks 64                #Number of processes
#SBATCH --nodes 1                  #Number of nodes
#SBATCH -t 300                 #Runtime in minutes
#SBATCH -p general   	      #Partition to submit to

#SBATCH --mem-per-cpu=200     #Memory per cpu in MB (see also --mem)
#SBATCH -o newfullOneNode256k.out     	      #File to which standard out will be written
#SBATCH -e newfullOneNode256k.err      	      #File to which standard err will be written
 
#
# Use modules to setup the runtime environment
#
. /etc/profile

module load centos6/openmpi-1.6.5_gcc-4.8.0 

make clean

make symmetric XFLAGS='-DP2P_NUM_THREADS=0 -DP2P_DECAY_ITERATOR=0'
#
# Execute the run
#

mpirun -np 1  ./symmetric 256000  -c 1 #-nocheck

mpirun -np 2  ./symmetric 256000  -c 1 #-nocheck

mpirun -np 4  ./symmetric 256000  -c 1 #-nocheck

mpirun -np 4  ./symmetric 256000  -c 2 #-nocheck

mpirun -np 8  ./symmetric 256000  -c 1 #-nocheck

mpirun -np 8  ./symmetric 256000  -c 2 #-nocheck

mpirun -np 16  ./symmetric 256000  -c 1 #-nocheck

mpirun -np 16  ./symmetric 256000  -c 2 #-nocheck

mpirun -np 16  ./symmetric 256000  -c 4 #-nocheck

mpirun -np 32  ./symmetric 256000  -c 1 #-nocheck

mpirun -np 32  ./symmetric 256000  -c 2 #-nocheck

mpirun -np 32  ./symmetric 256000  -c 4 #-nocheck

mpirun -np 64  ./symmetric 256000  -c 1 #-nocheck

mpirun -np 64  ./symmetric 256000  -c 2 #-nocheck

mpirun -np 64  ./symmetric 256000  -c 4 #-nocheck

mpirun -np 64  ./symmetric 256000  -c 8 #-nocheck