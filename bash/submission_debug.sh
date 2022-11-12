#!/bin/bash -l

#  R MPI parallel job

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=0:40:0

# Request 1 gigabyte of RAM per process.
#$ -l mem=150M

# Select the MPI parallel environment with 32 processes
#$ -pe mpi 50

# Load R/GRASS environment
echo "running init.sh script..."
source /home/tcrnbgh/Scratch/visual_prominence/bash/init.sh
echo "done!"

# Set the working directory to somewhere in your scratch space.  This is
# necessary because the compute nodes cannot write to your $HOME
# NOTE: this directory must exist.
# set working dir
#$ -wd /home/tcrnbgh/Scratch/visual_prominence

# Run our MPI job. GERun is our wrapper for mpirun, which launches MPI jobs  
echo "running gerun..."
gerun /home/tcrnbgh/RMPISNOW_bgh < /home/tcrnbgh/Scratch/visual_prominence/rscript/analysis_debug.R
echo "done!"