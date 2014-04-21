#!/bin/bash

make clean
make profile_p2p XFLAGS='-DP2P_BLOCK_SIZE=1'
echo "BLOCK_SIZE = 1"
./profile_p2p

make clean
make profile_p2p XFLAGS='-DP2P_BLOCK_SIZE=2'
echo "BLOCK_SIZE = 2"
./profile_p2p

make clean
make profile_p2p XFLAGS='-DP2P_BLOCK_SIZE=4'
echo "BLOCK_SIZE = 4"
./profile_p2p

make clean
make profile_p2p XFLAGS='-DP2P_BLOCK_SIZE=8'
echo "BLOCK_SIZE = 8"
./profile_p2p

make clean
make profile_p2p XFLAGS='-DP2P_BLOCK_SIZE=16'
echo "BLOCK_SIZE = 16"
./profile_p2p

make clean
make profile_p2p XFLAGS='-DP2P_BLOCK_SIZE=32'
echo "BLOCK_SIZE = 32"
./profile_p2p

make clean
make profile_p2p XFLAGS='-DP2P_BLOCK_SIZE=64'
echo "BLOCK_SIZE = 64"
./profile_p2p

make clean
make profile_p2p XFLAGS='-DP2P_BLOCK_SIZE=128'
echo "BLOCK_SIZE = 128"
./profile_p2p

make clean
make profile_p2p XFLAGS='-DP2P_BLOCK_SIZE=256'
echo "BLOCK_SIZE = 256"
./profile_p2p

make clean
make profile_p2p XFLAGS='-DP2P_BLOCK_SIZE=512'
echo "BLOCK_SIZE = 512"
./profile_p2p

make clean
make profile_p2p XFLAGS='-DP2P_BLOCK_SIZE=1024'
echo "BLOCK_SIZE = 1024"
./profile_p2p
