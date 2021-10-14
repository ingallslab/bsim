#!/bin/bash

# Shell script for setting up BSim environment on Graham compute cluster.
# Newest available version of java (13.0.2) is loaded in place of default (1.8.0).
# Creates alias for building bsim by executing 'bsim-build' in command line.

echo "Setting up BSim environment..."
module load ant
module load openmpi #change this to most updated environment
module load mpi4py #needed for multiprocessing in python and AstroABC code
module load StdEnv/2020 
module load scipy-stack/2020b
module load opencv
#module load java/13.0.2 #loaded AFTER StdEnv which make 13.0.2 available over 13.0.1 


export BSIM_DIR=$PWD
alias bsim-build='ant -f $BSIM_DIR/bsim-build-tree.xml'

# Load the python libraries and packages needed by ABC computations
# Should be stored in /path/to/bsim-hpc-package/python-env.sh 
if [[ ${HPC_BSIM+x} ]]; then
  source $HPC_BSIM/python-env.sh
else
  echo -e "Please install bsim-hpc-package (https://github.com/ingallslab/bsim-hpc-package) if you want to perform ABC parameter inference with BSim.\n"

echo -e "Done.\n"
