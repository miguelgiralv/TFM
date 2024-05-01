plink --freq --bfile <input> --out <freq-file>
# seguimos los pasos de https://imputationserver.readthedocs.io/en/latest/prepare-your-data/ 
wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

# ejecutamos plink en wsl:
"/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink" --freq --bfile \
"/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_final" \
--out "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/imputation/freq-file"
# sustituir por nombres de rutas bien más adelante para ejecutar las rutas
cp "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_final.bim" /home/miguel/impute
cp "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/imputation/freq-file"* /home/miguel/impute

perl HRC-1000G-check-bim.pl -b "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_final.bim" \
-f "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/imputation/freq-file" \
-r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h sh Run-plink.sh

# convertir a VCF:
plink --bfile [filename prefix] --recode vcf --out [VCF prefix]
"/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink" --bfile \
"/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_final" \
--recode vcf --out "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/imputation/GSE33528"

# bcf tools sort:
bcftools sort "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/imputation/GSE33528.vcf" \
-Oz -o /home/miguel/impute/study.vcf.gz

# ponemos el script de python:
python2 checkVCF.py -r human_g1k_v37.fasta -o out mystudy_chr1.vcf.gz

VcfCooker
