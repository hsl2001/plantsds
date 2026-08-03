#!/usr/bin/env bash
set -e

echo "=========================================================="
echo " Preparing SEDEF in WSL"
echo "=========================================================="

if [ -d "sedef" ]; then
    echo "Directory 'sedef' already exists. Pulling latest changes..."
    cd sedef
    git pull
else
    echo "Cloning SEDEF repository..."
    git clone https://github.com/vpc-ccg/sedef.git
    cd sedef
fi

echo "Building SEDEF..."
make

echo "=========================================================="
echo " SEDEF is ready at $(pwd)/sedef"
echo "=========================================================="
