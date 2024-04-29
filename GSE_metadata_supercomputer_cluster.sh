#!/bin/bash

#SBATCH -J analisis_metadata_gse_miguel
#SBATCH -p day
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o ./output_jobs/slurm.%N.%j.out
#SBATCH -e ./output_jobs/slurm.%N.%j.err

module load R

# Set local R library directory
export R_LIBS_USER=~/R/library

# Install BiocManager and GEOquery packages
R -e 'install.packages("BiocManager", repos="https://cran.r-project.org")'
R -e 'BiocManager::install(c("GEOquery"))'

# Run R script
#Rscript /home/mgiralv/tfm/gse/GSE_metadata_cluster.R

echo "done!"

