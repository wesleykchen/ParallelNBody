#!/bin/bash
#
#$ -cwd
##$ -j y
#$ -S /bin/bash
#$ -M wesleychen@seas.harvard.edu
#$ -pe orte 8
#$ -o output/clusterScript.out
#$ -e output/clusterScript.err
#
 
#
# Use modules to setup the runtime environment
#
. /etc/profile
module load compilers/intel/11.1
module load mpi/openmpi/1.4.2/intel
module load courses/cs205/2013
 
#
# Execute the run
#

mpirun -np $NSLOTS ./symmetric 1000 -c 2
