
path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"


# los cromosomas imputados los convertimos a plink
for chr in {1..22} X Y; do
    $path/software/plink --vcf $path/results/imputado/extracted/chr${chr}.dose.vcf.gz --make-bed --out $path/results/imputado/extracted/plink/imputed_data_chr${chr}
done


# calculamos los PRS (### prueba con genoma sin imputar ###) para PGS000054
$path/software/plink --bfile $path/results/plink_data/binary/processed/GSE33528_qc_cr_in \
--score $path/data/PRS/PGS000054_hmPOS_GRCh37.txt 1 4 5 header --out $path/results/PRS/PGS000054_PRS_results.txt

# para PGS004146
$path/software/plink --bfile $path/results/plink_data/binary/processed/GSE33528_qc_cr_in \
--score $path/data/PRS/PGS004146_hmPOS_GRCh37_2.txt 1 4 6 header --out $path/results/PRS/PGS004146_hmPOS_GRCh37

# Preparar Covariates File (covariates.profile) en r:

# ejecutar la regresion logistica
$path/software/plink --bfile $path/results/plink_data/binary/processed/GSE33528_qc_cr_in \
--covar $path/results/PRS/PGS004146_PRS_results.profile \
--covar-name SCORE \
--logistic \
--out $path/results/PRS/logistic_regression_with_PRS


