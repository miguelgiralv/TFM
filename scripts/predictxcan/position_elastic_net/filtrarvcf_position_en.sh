##########nuevo planteamiento 2:
path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

# 1. liftover de predictdb 
# hemos creado un archivo que contiene todas las posiciones para liftover
$path/software/liftOver \
    $path/results/imputado/extracted/all_positions_mashr.txt \
    $path/data/predictXcan/hg38ToHg19.over.chain.gz \
    $path/results/imputado/extracted/all_positions_mashr_hg19.txt \
    $path/results/imputado/extracted/all_positions_mashr_hg19_unlifted.txt

# 2. filtrado con las posiciones de predictdb
bcftools view \
-R $path/results/imputado/extracted/all_positions_mashr_hg19_nochr.txt \
$path/results/imputado/extracted/chr10.dose.vcf.gz \
-o $path/results/imputado/extracted/filtered/chr10.dose_filter_mashr.vcf.gz

# ahora filtramos todos los cromosomas: 
for chr in {2..5} ; do
bcftools view \
-R $path/results/imputado/extracted/all_positions_mashr_hg19_nochr.txt \
$path/results/imputado/extracted/chr${chr}.dose.vcf.gz \
-o $path/results/imputado/extracted/filtered/mashr/chr${chr}.dose_filter_mashr.vcf.gz
echo "${chr} filtrado!"
done

for chr in {6..10} ; do
bcftools view \
-R $path/results/imputado/extracted/all_positions_mashr_hg19_nochr.txt \
$path/results/imputado/extracted/chr${chr}.dose.vcf.gz \
-o $path/results/imputado/extracted/filtered/mashr/chr${chr}.dose_filter_mashr.vcf.gz
echo "${chr} filtrado!"
done

for chr in {10..15} ; do
bcftools view \
-R $path/results/imputado/extracted/all_positions_mashr_hg19_nochr.txt \
$path/results/imputado/extracted/chr${chr}.dose.vcf.gz \
-o $path/results/imputado/extracted/filtered/mashr/chr${chr}.dose_filter_mashr.vcf.gz
echo "${chr} filtrado!"
done

for chr in {16..20} ; do
bcftools view \
-R $path/results/imputado/extracted/all_positions_mashr_hg19_nochr.txt \
$path/results/imputado/extracted/chr${chr}.dose.vcf.gz \
-o $path/results/imputado/extracted/filtered/mashr/chr${chr}.dose_filter_mashr.vcf.gz
echo "${chr} filtrado!"
done

for chr in {21..22} X; do
bcftools view \
-R $path/results/imputado/extracted/all_positions_mashr_hg19_nochr.txt \
$path/results/imputado/extracted/chr${chr}.dose.vcf.gz \
-o $path/results/imputado/extracted/filtered/mashr/chr${chr}.dose_filter_mashr.vcf.gz
echo "${chr} filtrado!"
done

for chr in {2..5} {6..10} {10..15} {16..20} {21..22} X; do
bcftools view \
-R $path/results/imputado/extracted/all_positions_mashr_hg19_nochr.txt \
$path/results/imputado/extracted/chr5.dose.vcf.gz \
-o $path/results/imputado/extracted/filtered/mashr/chr5.dose_filter_mashr.vcf.gz
echo "${chr} filtrado!"
done


############# Concatenamos todos los cromosomas:
# por si hay que indexar antes:
for chr in {1..22} X; do
    bcftools index $path/results/imputado/extracted/filtered/mashr/chr${chr}.dose_filter_mashr.vcf.gz
    echo "${chr} indexado!"
done

bcftools concat -Oz -o $path/results/imputado/extracted/filtered/mashr/fullgenome_mashr.vcf.gz \
    $path/results/imputado/extracted/filtered/mashr/chr{1..22}.dose_filter_mashr.vcf.gz \
    $path/results/imputado/extracted/filtered/mashr/chrX.dose_filter_mashr.vcf.gz \

####### Finalmente aplicamos predictXcan con sangre
python3 $METAXCAN/Predict.py \
--model_db_path $DATA/models/gtex_v8_mashr/mashr_Whole_Blood.db \
--model_db_snp_key varID \
--vcf_genotypes $path/results/imputado/extracted/filtered/mashr/fullgenome_mashr.vcf.gz \
--vcf_mode imputed \
--liftover $DATA/hg19ToHg38.over.chain.gz \
--on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
--prediction_output $RESULTS/vcf_1000G_hg37_mashr/test_whole_Whole_Blood_7_predict.txt \
--prediction_summary_output $RESULTS/vcf_1000G_hg37_mashr/test_whole_Whole_Blood_7_summary.txt \
--verbosity 9 \
--throw


## aplicamos predictXcan en todos los tejidos
path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"
METAXCAN="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/software/MetaXcan"
DATA="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/data/predictXcan"
RESULTS="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/predictXcan/test"

{
for tissue_file in $DATA/predictXcan/mashr_models/*.db; do
    tissue=$(basename "$tissue_file" .db)
    python3 $METAXCAN/Predict.py \
    --model_db_path "$tissue_file" \
    --model_db_snp_key varID \
    --vcf_genotypes $path/results/imputado/extracted/filtered/mashr/fullgenome_mashr.vcf.gz \
    --vcf_mode imputed \
    --liftover $DATA/hg19ToHg38.over.chain.gz \
    --on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
    --prediction_output $RESULTS/vcf_1000G_hg37_mashr/${tissue}_predict.txt \
    --prediction_summary_output $RESULTS/vcf_1000G_hg37_mashr/${tissue}_summary.txt \
    --verbosity 9 \
    --throw 
done
} 2>&1 | tee "$RESULTS/vcf_1000G_hg37_mashr/log.txt"


## FUNCIONA LOOP POR TODOS LOS TEJIDOS:
path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"
METAXCAN="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/software/MetaXcan"
DATA="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/data/predictXcan"
RESULTS="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/predictXcan/test"

{
for tissue_file in "$DATA/mashr_models"/*.db; do
    echo "Processing file: $tissue_file"  # Debugging output
    tissue=$(basename "$tissue_file" .db)
    python3 "$METAXCAN/Predict.py" \
    --model_db_path "$tissue_file" \
    --model_db_snp_key varID \
    --vcf_genotypes "$path/results/imputado/extracted/filtered/mashr/fullgenome_mashr.vcf.gz" \
    --vcf_mode imputed \
    --liftover "$DATA/hg19ToHg38.over.chain.gz" \
    --on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
    --prediction_output "$RESULTS/vcf_1000G_hg37_mashr/${tissue}_predict.txt" \
    --prediction_summary_output "$RESULTS/vcf_1000G_hg37_mashr/${tissue}_summary.txt" \
    --verbosity 9 \
    --throw 
done
} 2>&1 | tee "$RESULTS/vcf_1000G_hg37_mashr/log.txt"

# LOOP POR PARTES
{
for tissue_file in "$DATA/mashr_models/1"/*.db; do
    echo "Processing file: $tissue_file"  # Debugging output
    tissue=$(basename "$tissue_file" .db)
    python3 "$METAXCAN/Predict.py" \
    --model_db_path "$tissue_file" \
    --model_db_snp_key varID \
    --vcf_genotypes "$path/results/imputado/extracted/filtered/mashr/fullgenome_mashr.vcf.gz" \
    --vcf_mode imputed \
    --liftover "$DATA/hg19ToHg38.over.chain.gz" \
    --on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
    --prediction_output "$RESULTS/vcf_1000G_hg37_mashr/${tissue}_predict.txt" \
    --prediction_summary_output "$RESULTS/vcf_1000G_hg37_mashr/${tissue}_summary.txt" \
    --verbosity 9 \
    --throw 
done
} 2>&1 | tee "$RESULTS/vcf_1000G_hg37_mashr/1/log.txt"


{
for tissue_file in "$DATA/mashr_models/5"/*.db; do
    echo "Processing file: $tissue_file"  
    tissue=$(basename "$tissue_file" .db)
    python3 "$METAXCAN/Predict.py" \
    --model_db_path "$tissue_file" \
    --model_db_snp_key varID \
    --vcf_genotypes "$path/results/imputado/extracted/filtered/mashr/fullgenome_mashr.vcf.gz" \
    --vcf_mode imputed \
    --liftover "$DATA/hg19ToHg38.over.chain.gz" \
    --on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
    --prediction_output "$RESULTS/vcf_1000G_hg37_mashr/${tissue}_predict.txt" \
    --prediction_summary_output "$RESULTS/vcf_1000G_hg37_mashr/${tissue}_summary.txt" \
    --verbosity 9 \
    --throw 
done
} 2>&1 | tee "$RESULTS/vcf_1000G_hg37_mashr/5/log.txt"
