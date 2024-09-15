#!/bin/bash
#SBATCH -J imputacion_vcf_miguel
#SBATCH -p day
#SBATCH -N 2
#SBATCH -n 4
#SBATCH -o ./output_jobs/slurm.%N.%j.out
#SBATCH -e ./output_jobs/slurm.%N.%j.err

# echo en el supercomputador: salloc -J imputacion_miguel -p day -N 1 -c 4
perl HRC-1000G-check-bim-NoReadKey.pl -b "/home/mgiralv/tfm/imputation/GSE33528_qc_cr_in.bim" \
-f "/home/mgiralv/tfm/imputation/freq-file" \
-r "/home/mgiralv/tfm/imputation/freq-file/1000GP_Phase3_combined.legend" -g \
sh Run-plink.sh \

echo "done!"
