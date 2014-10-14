#!/bin/bash
#BATCH -n 64                #Number of cores
#SBATCH -t 10                 #Runtime in minutes
#SBATCH -p general   	      #Partition to submit to

#SBATCH --mem-per-cpu=2000     #Memory per cpu in MB (see also --mem)
#SBATCH -o t1_64_4.out     	      #File to which standard out will be written
#SBATCH -e t1_64_4.err      	      #File to which standard err will be written
 
#
# Use modules to setup the runtime environment
#
. /etc/profile
#module load centos6/openmpi-1.7.2_intel-13.0.079
#module load centos6/fftw-3.3.3_openmpi-1.6.4_gcc-4.8.0
module load centos6/openmpi-1.6.5_gcc-4.8.0 
#
# Execute the run
#

mpirun -np 32  ./symmetric 25600  -c 4

mpirun -np 64  ./symmetric 25600  -c 4

mpirun -np 32  ./teamscatter 25600  -c 4

mpirun -np 64  ./teamscatter 25600  -c 4

#mpirun -np 64 ./symmetric 25600 -c 2 -nocheck

#mpirun -np 128 ./symmetric 25600 -c 4 -nocheck

#mpirun -np 128 ./symmetric 25600 -c 8 -nocheck

#mpirun -np 128 ./teamscatter 25600 -c 1 -nocheck

#mpirun -np 128 ./teamscatter 25600 -c 2 -nocheck

#mpirun -np 128 ./teamscatter 25600 -c 4 -nocheck

#mpirun -np 128 ./teamscatter 25600 -c 8 -nocheck