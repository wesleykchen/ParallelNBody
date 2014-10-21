#!/bin/bash

OUTA=data/p2p_A.out
OUTO=data/p2p_O.out
OUTD=data/p2p_D.out

for B in 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576 2097152 4194304
do
  make clean
  make profile_p2p XFLAGS="-DP2P_BLOCK_SIZE=$B -DP2P_NUM_THREADS=0"
  ./profile_p2p 'A' 2>&1 | tee -a ${OUTA}
	echo -e "\n" 2>&1 | tee -a ${OUTA}
  ./profile_p2p 'O' 2>&1 | tee -a ${OUTO}
	echo -e "\n" 2>&1 | tee -a ${OUTO}
  ./profile_p2p 'D' 2>&1 | tee -a ${OUTD}
	echo -e "\n" 2>&1 | tee -a ${OUTD}
done
