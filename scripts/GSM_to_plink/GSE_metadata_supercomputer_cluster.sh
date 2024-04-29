#!/bin/bash

#SBATCH -J analisis_metadata_gse_miguel
#SBATCH -p day
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o ./output_jobs/slurm.%N.%j.out
#SBATCH -e ./output_jobs/slurm.%N.%j.err

module load R

# Directorio de R a usar para las librerias 
export R_LIBS_USER=~/R/library

# Instalar BiocManager y GEOquery packages
#R -e 'install.packages("BiocManager", repos="https://cran.r-project.org")'
#R -e 'BiocManager::install(c("GEOquery"))'

# correr script
#Rscript /home/mgiralv/tfm/gse/GSE_metadata_cluster.R

echo "done!"

