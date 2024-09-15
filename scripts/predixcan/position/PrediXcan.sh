path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

# definimos nuevas rutas
path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"
METAXCAN="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/software/MetaXcan"
DATA="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/data/predictXcan"
RESULTS="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/predictXcan/test"

## aplicamos predictXcan en todos los tejidos

# completo
{
for tissue_file in "$DATA/mashr_models"/*.db; do
    echo "Procesando: $tissue_file"  
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


# Alternativa de computacion paralela
# Loop por partes para ejecutar en varias consolas de wsl en paralelo: (dividimos los tejidos en 5 carpetas)
{
for tissue_file in "$DATA/mashr_models/1"/*.db; do
    echo "Procesando: $tissue_file"  
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
for tissue_file in "$DATA/mashr_models/2"/*.db; do
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
} 2>&1 | tee "$RESULTS/vcf_1000G_hg37_mashr/2/log.txt"

{
for tissue_file in "$DATA/mashr_models/3"/*.db; do
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
} 2>&1 | tee "$RESULTS/vcf_1000G_hg37_mashr/3/log.txt"

{
for tissue_file in "$DATA/mashr_models/4"/*.db; do
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
} 2>&1 | tee "$RESULTS/vcf_1000G_hg37_mashr/4/log.txt"


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
