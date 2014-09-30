#!/bin/bash
#SBATCH -n 4                 						#Number of cores
#SBATCH -t 5                 						#Runtime in minutes
#SBATCH -p general   								#Partition to submit to
#SBATCH --mem-per-cpu=100   						#Memory per cpu in MB (see also --mem)
#SBATCH -o clusterScript.out     					#File to which standard out will be written
#SBATCH -e clusterScript.err      					#File to which standard err will be written
 
#
# Use modules to setup the runtime environment
#
. /etc/profile
module load hpc/intel-compilers-13.0.079 
module load hpc/openmpi-1.6.2_gcc-4.7.2
module load centos6/openmpi-1.7.2_intel-13.0.079
 
#
# Execute the run
#

mpirun -np 4 ./symmetric 1000 -c 2

mpirun -np 2 ./symmetric 1000 -c 1

mpirun -np 1 ./symmetric 1000 -c 1
