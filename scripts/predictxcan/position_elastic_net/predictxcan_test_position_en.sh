path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"


METAXCAN="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/software/MetaXcan"
DATA="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/data/predictXcan"
RESULTS="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/predictXcan/test"

# solo un tejido:
python3 $METAXCAN/Predict.py \
--model_db_path $DATA/elastic_net_models/en_Adipose_Subcutaneous.db \
--model_db_snp_key varID \
--vcf_genotypes $path/results/imputado/extracted/filtered/elastic_net/chr1.dose_filter_en.vcf.gz \
--vcf_mode imputed \
--liftover $DATA/hg19ToHg38.over.chain.gz \
--on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
--prediction_output $RESULTS/vcf_1000G_hg37_en/Adipose_Subcutaneous.txt \
--prediction_summary_output $RESULTS/vcf_1000G_hg37_en/Adipose_Subcutaneous_summary.txt \
--verbosity 9 \
--throw


#solo un tejido ()
python3 $METAXCAN/Predict.py \
--model_db_path $DATA/models/gtex_v8_en/mashr_Whole_Blood.db \
--model_db_snp_key varID \
--vcf_genotypes $path/results/imputado/extracted/filtered/elastic_net/chr1.dose_filter_en.vcf.gz \
--vcf_mode imputed \
--liftover $DATA/hg19ToHg38.over.chain.gz \
--on_the_fly_mapping METADATA "{}_{}_{}_{}"
--prediction_output $RESULTS/vcf_1000G_hg37_mashr/test_whole_Whole_Blood_7_predict.txt \
--prediction_summary_output $RESULTS/vcf_1000G_hg37_mashr/test_whole_Whole_Blood_7_summary.txt \
--verbosity 9 \
--throw


# todos los tejidos:
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





# ejemplo 1: VCF hg19-based genotype, on GTEx v8 Elastic Net Model (funciona con los cromosomas filtrados y en rsid)
printf "Predict expression\n\n"
python3 $METAXCAN/Predict.py \
--model_db_path $DATA/models/gtex_v8_en/en_Whole_Blood.db \
--vcf_genotypes $path/results/imputado/extracted/filtered/chr1.dose_rsid_filter.vcf.gz  \
--vcf_mode genotyped \
--prediction_output $RESULTS/vcf_1000G_hg37_en/test_Whole_Blood__predict.txt \
--prediction_summary_output $RESULTS/vcf_1000G_hg37_en/test_Whole_Blood__summary.txt \
--verbosity 9 \
--throw

# prediccion con posiciones en hg19 (sin rsid) haciendo liftover

python3 $METAXCAN/Predict.py \
--model_db_path $DATA/models/gtex_v8_mashr/mashr_Whole_Blood.db \
--model_db_snp_key varID \
--vcf_genotypes $path/results/imputado/extracted/chr1.dose.vcf.gz \
--vcf_mode genotyped \
--liftover $DATA/hg19ToHg38.over.chain.gz \
--on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
--prediction_output $RESULTS/vcf_1000G_hg37_mashr/test2_Whole_Blood__predict.txt \
--prediction_summary_output $RESULTS/vcf_1000G_hg37_mashr/test2_Whole_Blood__summary.txt \
--verbosity 9 \
--throw

## 3 opciones:

#########1. filtrando predictdb sin rsid

# 1. liftover de predictdb 

# 2. filtrado con las posiciones

# 3. predicción con posiciones


#########2. dividiendo por muestras sin filtrar

# 