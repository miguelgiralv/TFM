path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"


METAXCAN="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/software/MetaXcan/"
DATA="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/data/predictXcan"
RESULTS="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/predictXcan/test"

#probar con hg19 filtrado con las posiciones pero sin rsid









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