#!/bin/bash

make clean --silent
echo "BLOCK_SIZE = 1"
make profile_p2p CFLAGS='-DP2P_BLOCK_SIZE=1' --silent
./profile_p2p

make clean --silent
echo "BLOCK_SIZE = 2"
make profile_p2p CFLAGS='-DP2P_BLOCK_SIZE=2' --silent
./profile_p2p

make clean --silent
echo "BLOCK_SIZE = 4"
make profile_p2p CFLAGS='-DP2P_BLOCK_SIZE=4' --silent
./profile_p2p

make clean --silent
echo "BLOCK_SIZE = 8"
make profile_p2p CFLAGS='-DP2P_BLOCK_SIZE=8' --silent
./profile_p2p

make clean --silent
echo "BLOCK_SIZE = 16"
make profile_p2p CFLAGS='-DP2P_BLOCK_SIZE=16' --silent
./profile_p2p

make clean --silent
echo "BLOCK_SIZE = 32"
make profile_p2p CFLAGS='-DP2P_BLOCK_SIZE=32' --silent
./profile_p2p

make clean --silent
echo "BLOCK_SIZE = 64"
make profile_p2p CFLAGS='-DP2P_BLOCK_SIZE=64' --silent
./profile_p2p

make clean --silent
echo "BLOCK_SIZE = 128"
make profile_p2p CFLAGS='-DP2P_BLOCK_SIZE=128' --silent
./profile_p2p

make clean --silent
echo "BLOCK_SIZE = 256"
make profile_p2p CFLAGS='-DP2P_BLOCK_SIZE=256' --silent
./profile_p2p

make clean --silent
echo "BLOCK_SIZE = 512"
make profile_p2p CFLAGS='-DP2P_BLOCK_SIZE=512' --silent
./profile_p2p

make clean --silent
echo "BLOCK_SIZE = 1024"
make profile_p2p CFLAGS='-DP2P_BLOCK_SIZE=1024' --silent
./profile_p2p
