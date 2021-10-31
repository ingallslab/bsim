#!/bin/bash

# Shell script for setting up BSim environment on Graham compute cluster.
# Newest available version of java (13.0.2) is loaded in place of default (1.8.0).
# Creates alias for building bsim by executing 'bsim-build' in command line.

echo "Setting up BSim environment..."
module load ant
module load openmpi
module load mpi4py 
module load StdEnv/2020 
module load scipy-stack
module load opencv



# Load the python libraries and packages needed by ABC computations
# Should be stored in /path/to/bsim-hpc-package/python-env.sh 
if [[ ${HPC_BSIM+x} ]]; then
  export BSIM_DIR=$HPC_BSIM/bsim
  alias bsim-build='ant -f $BSIM_DIR/bsim-build-tree.xml'
  source $HPC_BSIM/python-env.sh
else
  echo -e "Please install bsim-hpc-package (https://github.com/ingallslab/bsim-hpc-package) if you want to perform ABC parameter inference with BSim or use high performance computing resources on ComputeCanada infrastructure..\n"
fi

echo -e "Done.\n"
