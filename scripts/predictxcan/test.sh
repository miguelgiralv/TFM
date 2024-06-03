path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/test_master"

# descargamos predictXcan y adecuamos el entorno de python siguiendo este tutorial:
# https://github.com/hakyimlab/MetaXcan 

# primero concatenamos los vcf de los cromosomas que hemos obtenido del servidor de michigan para tener un único vcf:

# los indexamos primero
path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

for chr in {1..22} X; do
    bcftools index $path/results/imputado/extracted/chr${chr}.dose.vcf.gz
    echo "${chr} indexado!"
done

bcftools concat -Oz -o $path/results/imputado/extracted/fullgenomesdose.vcf.gz \
    $path/results/imputado/extracted/chr{1..22}.dose.vcf.gz \
    $path/results/imputado/extracted/chrX.dose.vcf.gz



METAXCAN="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/software/MetaXcan/"
DATA="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/data/predictXcan"
RESULTS="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/predictXcan/test"

# ejemplo
printf "Predict expression\n\n"
python3 $METAXCAN/Predict.py \
--model_db_path $DATA/models/gtex_v8_en/en_Whole_Blood.db \
--vcf_genotypes $DATA/1000G_hg37/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--vcf_mode genotyped \
--prediction_output $RESULTS/vcf_1000G_hg37_en/Whole_Blood__predict.txt \
--prediction_summary_output $RESULTS/vcf_1000G_hg37_en/Whole_Blood__summary.txt \
--verbosity 9 \
--throw


# test con mi vcf:
printf "Predict expression\n\n"
python3 $METAXCAN/Predict.py \
--model_db_path $DATA/models/gtex_v8_en/en_Whole_Blood.db \
--model_db_snp_key varID \
--vcf_genotypes $path/results/imputado/extracted/fullgenomesdose.vcf.gz \
--vcf_mode imputed \
--prediction_output $RESULTS/test_Whole_Blood__predict.txt \
--prediction_summary_output $RESULTS/test_Whole_Blood__summary.txt \
--verbosity 9 \
--throw



# ver la informacion de los archivos .info

# Assuming your VCF file is named input.vcf and SNP list is snp_list.txt

# Sort and index the SNP list file (if not already sorted and indexed)
sort -u $path/results/imputado/extracted/snps_predictdb.txt > $path/results/imputado/extracted/snps_predictdb_sorted.txt

awk -F':' '{print $1"\t"$2}' $path/results/imputado/extracted/snps_predictdb.txt > $path/results/imputado/extracted/snps_predictdb_corrected.txt

awk -F: '{print $1"\t"$2}' $path/results/imputado/extracted/snps_predictdb_sorted.txt > $path/results/imputado/extracted/snps_predictdb_corrected.txt

awk -F: '{print $1"\t"$2}' snps.txt > snps_processed.txt

bcftools index $path/results/imputado/extracted/fullgenomesdose.vcf.gz

sed 's/:/\t/g' $path/results/imputado/extracted/snps_predictdb.txt > $path/results/imputado/extracted/snps_predictdb_corrected.txt



# Filter the VCF file based on the SNP list
bcftools view -T $path/results/imputado/extracted/snps_predictdb_sorted.txt  $path/results/imputado/extracted/fullgenomesdose.vcf.gz -o output_preditdb.vcf

# vamos a probar a filtrar cromosoma a cromosoma y despues unirlo todo en un unico vcf
bcftools view -r 10:295356 $path/results/imputado/extracted/chr1.dose.vcf.gz -o $path/results/imputado/extracted/chr1.dose_preditdb.vcf.gz

# separado de tabuladores
bcftools view -R $path/results/imputado/extracted/snps_predictdb_corrected.txt $path/results/imputado/extracted/chr1.dose.vcf.gz -o $path/results/imputado/extracted/chr1.dose_preditdb.vcf.gz

bcftools view -t  $path/results/imputado/extracted/chr1.dose.vcf.gz -o $path/results/imputado/extracted/chr1.dose_preditdb.vcf.gz
