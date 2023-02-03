#!/usr/bin/env bash
#SBATCH -A C3SE2023-2-2 # Project name
#SBATCH -J fd_LDA # Name of the job
#SBATCH -N 1 # Amount of nodes
#SBATCH -n 1 # Amount of cores per node
#SBATCH -t 1:00:00 # Maximum run time
#SBATCH -o fd_LDA.out # All non-error output
#SBATCH -e fd_LDA_err.out # All error-related output

# Unload modules
module purge

# Load necessary modules
module load ASE/3.22.1-foss-2022a
module load GPAW/22.8.0-foss-2022a

# Run the program
srun gpaw python fd_LDA.py
