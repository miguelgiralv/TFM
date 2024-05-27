#!/bin/bash
#SBATCH -J vcf_filtrarmiguel
#SBATCH -p day
#SBATCH -N 2
#SBATCH -n 4
#SBATCH -o ./output_jobs/slurm.%N.%j.out
#SBATCH -e ./output_jobs/slurm.%N.%j.err

module load miniconda3
conda activate bcftools_env


bcftools index /home/mgiralv/tfm/predictxcan/fullgenomesdose.vcf.gz

echo "indexado!"

# Filter the VCF file based on the SNP list
bcftools view -T /home/mgiralv/tfm/predictxcan/snps_predictdb_sorted.txt  /home/mgiralv/tfm/predictxcan/fullgenomesdose.vcf.gz -o /home/mgiralv/tfm/predictxcan/output_preditdb.vcf


echo "done!"
