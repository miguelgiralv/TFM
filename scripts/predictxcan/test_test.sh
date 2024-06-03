#!/bin/bash
#SBATCH -J vcf_filtrarmiguel
#SBATCH -p day
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o ./output_jobs/slurm.%N.%j.out
#SBATCH -e ./output_jobs/slurm.%N.%j.err

conda activate bcftools_env

echo "test"

