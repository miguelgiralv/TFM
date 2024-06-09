##########nuevo planteamiento 2:

# 1. liftover de predictdb 
# hemos creado un archivo que contiene todas las posiciones para liftover
$path/software/liftOver \
    $path/results/imputado/extracted/all_positions_mashr.txt \
    $path/data/predictXcan/hg38ToHg19.over.chain.gz \
    $path/results/imputado/extracted/all_positions_mashr_hg19.txt \
    $path/results/imputado/extracted/all_positions_mashr_hg19_unlifted.txt

# 2. filtrado con las posiciones de predictdb (FUNCIONA!))
bcftools view \
-R $path/results/imputado/extracted/all_positions_mashr_hg19_nochr.txt \
$path/results/imputado/extracted/chr1.dose.vcf.gz \
-o $path/results/imputado/extracted/filtered/chr1.dose_filter_mashr.vcf.gz

# ahora filtramos todos los cromosomas: 
for chr in {2..5} ; do
bcftools view \
-R $path/results/imputado/extracted/all_positions_mashr_hg19_nochr.txt \
$path/results/imputado/extracted/chr${chr}.dose.vcf.gz \
-o $path/results/imputado/extracted/filtered/mashr/chr${chr}.dose_filter_mashr.vcf.gz
echo "${chr} filtrado!"
done



# 3. predicción con posiciones filtradas (nos ahorramos lo de los rsid)
path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"
METAXCAN="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/software/MetaXcan/"
DATA="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/data/predictXcan"
RESULTS="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/predictXcan/test"

#probar con hg19 filtrado con las posiciones pero sin rsid

# prediccion con posiciones en hg19 (sin rsid) haciendo liftover
python3 $METAXCAN/Predict.py \
--model_db_path $DATA/models/gtex_v8_mashr/mashr_Whole_Blood.db \
--model_db_snp_key varID \
--vcf_genotypes $path/results/imputado/extracted/filtered/chr1.dose_filter_mashr.vcf.gz \
--vcf_mode genotyped \
--liftover $DATA/hg19ToHg38.over.chain.gz \
--on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
--prediction_output $RESULTS/vcf_1000G_hg37_mashr/test2_Whole_Blood__predict.txt \
--prediction_summary_output $RESULTS/vcf_1000G_hg37_mashr/test2_Whole_Blood__summary.txt \
--verbosity 9 \
--throw






















############### primero debemos anotar los vcf con sus rsid
# descargamos el vcf de referencia de hg 19: https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/ 
path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

# el vcf debemos procesarlo, ya que su columna viene en formato chr1:
sed 's/ID=chr/ID=/g' $path/results/imputation/1000genomas/header.vcf > $path/results/imputation/1000genomas/new_header.vcf
bcftools reheader -h $path/results/imputation/1000genomas/new_header.vcf -o $path/results/imputation/1000genomas/00-All_2.vcf.gz $path/results/imputation/1000genomas/00-All.vcf.gz
bcftools index $path/results/imputation/1000genomas/00-All_2.vcf.gz
bcftools annotate --rename-chrs $path/results/imputation/1000genomas/rename_chr.txt $path/results/imputation/1000genomas/00-All_2.vcf.gz -o $path/results/imputation/1000genomas/00-All_3.vcf.gz
bcftools index $path/results/imputation/1000genomas/00-All_3.vcf.gz

#una vez procesado podemos anotar los rsid con bcftools:

bcftools annotate \
    -a $path/results/imputation/1000genomas/00-All_3.vcf.gz \
    -c CHROM,POS,REF,ALT,ID \
    -o $path/results/imputado/extracted/rsid/chr2.dose_rsid.vcf.gz \
    $path/results/imputado/extracted/chr2.dose.vcf.gz

# vemos los resultados (funciona)
bcftools view $path/results/imputation/1000genomas/chr1.dose_rsid.vcf.gz | less

# ahora hacemos un bucle for para anotar todos los cromosomas. Como tarda mucho haremos 4 bucles
for chr in {2..5} ; do
    bcftools annotate \
    -a $path/results/imputation/1000genomas/00-All_3.vcf.gz \
    -c CHROM,POS,REF,ALT,ID \
    -o $path/results/imputado/extracted/rsid/chr${chr}.dose_rsid.vcf.gz \
    $path/results/imputado/extracted/chr${chr}.dose.vcf.gz
    echo "${chr} anotado!"
done

for chr in {6..11} ; do
    bcftools annotate \
    -a $path/results/imputation/1000genomas/00-All_3.vcf.gz \
    -c CHROM,POS,REF,ALT,ID \
    -o $path/results/imputado/extracted/rsid/chr${chr}.dose_rsid.vcf.gz \
    $path/results/imputado/extracted/chr${chr}.dose.vcf.gz
    echo "${chr} anotado!"
done

for chr in {12..17} ; do
    bcftools annotate \
    -a $path/results/imputation/1000genomas/00-All_3.vcf.gz \
    -c CHROM,POS,REF,ALT,ID \
    -o $path/results/imputado/extracted/rsid/chr${chr}.dose_rsid.vcf.gz \
    $path/results/imputado/extracted/chr${chr}.dose.vcf.gz
    echo "${chr} anotado!"
done

for chr in {17..22} X; do
    bcftools annotate \
    -a $path/results/imputation/1000genomas/00-All_3.vcf.gz \
    -c CHROM,POS,REF,ALT,ID \
    -o $path/results/imputado/extracted/rsid/chr${chr}.dose_rsid.vcf.gz \
    $path/results/imputado/extracted/chr${chr}.dose.vcf.gz
    echo "${chr} anotado!"
done

#uno a uno:
bcftools annotate \
-a $path/results/imputation/1000genomas/00-All_3.vcf.gz \
-c CHROM,POS,REF,ALT,ID \
-o $path/results/imputado/extracted/rsid/chr2.dose_rsid.vcf.gz \
$path/results/imputado/extracted/chr2.dose.vcf.gz
echo "${chr} anotado!"


############# ahora filtramos los rsid con las variantes de predictdb:

for chr in {1..5} ; do
    bcftools view \
    -R $path/data/predictXcan/elastic_net_models/all_rsids.txt \
    -o $path/results/imputado/extracted/rsid/chr${chr}.dose_filter_rsid.vcf.gz \
    $path/results/imputado/extracted/rsid/chr${chr}.dose_rsid.vcf.gz
    echo "${chr} filtrado!"
done

#FUNCIONA PARA FILTRAR, AUNQUE TARDA
bcftools view -i ID=@$path/data/predictXcan/elastic_net_models/all_rsids.txt \
$path/results/imputado/extracted/rsid/chr1.dose_rsid.vcf.gz \
-o $path/results/imputado/extracted/filtered/chr1.dose_rsid_filter.vcf.gz

bcftools view \
-R $path/data/predictXcan/elastic_net_models/all_rsids.txt \
$path/results/imputado/extracted/rsid/chr1.dose_rsid.vcf.gz \
-o $path/results/imputado/extracted/filtered/chr1.dose_filter_rsid.vcf.gz

#VEMOS EL NUMERO DE SNPS, que habia y luego tenemos:

for chr in {1..5} ; do
    n_${chr}_og = bcftools view -v snps sample.vcf | grep -v "^#" | wc -l 
    n_${chr}_rsid = bcftools view -v $path/results/imputado/extracted/rsid/chr1.dose_rsid.vcf.gz | grep -v "^#" | wc -l 
    n_${chr}_filt = bcftools view -v $path/results/imputado/extracted/filtered/chr1.dose_filter_rsid.vcf.gz | grep -v "^#" | wc -l 

    n_1_og = bcftools view -v snps sample.vcf | grep -v "^#" | wc -l 
    n_1_rsid = bcftools view -v $path/results/imputado/extracted/rsid/chr1.dose_rsid.vcf.gz | grep -v "^#" | wc -l 
    n_1_filt = bcftools view -v $path/results/imputado/extracted/filtered/chr1.dose_filter_rsid.vcf.gz | grep -v "^#" | wc -l 
    
    
bcftools view \
-R $path/data/predictXcan/elastic_net_models/all_rsids.txt \
-o $path/results/imputado/extracted/filtered/chr${chr}.dose_filter_rsid.vcf.gz \
$path/results/imputado/extracted/rsid/chr${chr}.dose_rsid.vcf.gz
echo "${chr} filtrado!"
done

################### prueba en el cromosoma 1. lo dividimos en dos partes y 

# lo partimos por la mitad de casos
bcftools query -l $path/results/imputado/extracted/chr1.dose.vcf.gz > $path/results/imputado/extracted/rsid/samples.txt
total_samples=$(wc -l < $path/results/imputado/extracted/rsid/samples.txt)
half_samples=$((total_samples / 2))
head -n $half_samples $path/results/imputado/extracted/rsid/samples.txt > $path/results/imputado/extracted/rsid/selected_samples.txt
tail -n +$((half_samples + 1)) $path/results/imputado/extracted/rsid/samples.txt > $path/results/imputado/extracted/rsid/remaining_samples.txt

bcftools view -S $path/results/imputado/extracted/rsid/selected_samples.txt -o $path/results/imputado/extracted/chr_samples/chr1_1.vcf.gz -O v $path/results/imputado/extracted/chr1.dose.vcf.gz
bcftools view -S $path/results/imputado/extracted/rsid/remaining_samples.txt -o $path/results/imputado/extracted/chr_samples/chr1_2.vcf.gz -O v $path/results/imputado/extracted/chr1.dose.vcf.gz

chr_samples

bcftools annotate \
-a $path/results/imputation/1000genomas/00-All_3.vcf.gz \
-c CHROM,POS,REF,ALT,ID \
-o $path/results/imputado/extracted/rsid/chr1_1.vcf.gz \
$path/results/imputado/extracted/chr_samples/chr1_1.vcf.gz
echo "${chr} anotado!"

bcftools view -R $path/data/predictXcan/elastic_net_models/all_rsids_chrom.txt $path/results/imputado/extracted/rsid/chr1.dose_rsid.vcf.gz -o $path/results/imputado/extracted/filtered/chr1.dose_rsid_filter.vcf.gz

################################PRUEBAS

bcftools view -i ID=@$path/data/predictXcan/elastic_net_models/all_rsids.txt $path/results/imputado/extracted/rsid/chr1.dose_rsid.vcf.gz \
-o $path/results/imputado/extracted/filtered/chr1.dose_rsid_filter.vcf.gz

bcftools view \
--samples-file $path/data/predictXcan/elastic_net_models/all_rsids.txt \
$path/results/imputado/extracted/rsid/chr1.dose_rsid.vcf.gz \
-Oz -o $path/results/imputado/extracted/filtered/chr1.dose_rsid_filter.vcf.gz


bcftools view -R $path/data/predictXcan/elastic_net_models/all_rsids.txt \
$path/results/imputado/extracted/rsid/chr1.dose_rsid.vcf.gz \
-Oz -o $path/results/imputado/extracted/filtered/chr1.dose_rsid_filter.vcf.gz


# visualizamos el nuevo 
bcftools view $path/results/imputado/extracted/filtered/chr1.dose_rsid_filter.vcf.gz | less

bcftools view -r 10:295356 $path/results/imputado/extracted/filtered/chr1.dose_rsid_filter.vcf.gz | less

bcftools view $DATA/1000G_hg37/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | less

############################# No revisado

#finalmente, juntamos los vcf procesados en un unico vcf que pasaramos por predictscan



#####################################

##########nuevo planteamiento 1

#1. liftover la bd de predictdb a hg19, (si no se puede filtrar por rsid)

#2 seleccionar con view solo esas posicioanes en todos los cromosomas

# 3 rsid anotar

# 4 unir todo 

#5predictdb

##########nuevo planteamiento 2:

# 1. liftover de predictdb 
# hemos creado un archivo para liftover
$path/software/liftOver \
    $path/results/imputado/extracted/snps_predictdb.bed \
    $path/data/predictXcan/hg38ToHg19.over.chain.gz \
    $path/results/imputado/extracted/snp_predictdb_hg19.bed \
    $path/results/imputado/extracted/snp_predictdb_hg19_unlifted.bed

$path/software/liftOver \
    -minMatch=0 \
    $path/results/imputado/extracted/snps_predict_chr.bed \
    $path/data/predictXcan/hg38ToHg19.over.chain.gz \
    $path/results/imputado/extracted/snp_predictdb_hg19.bed \
    $path/results/imputado/extracted/snp_predictdb_hg19_unlifted.bed


\data\predictXcan

# 2. filtrado con las posiciones de predictdb (NO FILTRA NADA!!!! :( ))
bcftools view \
-R $path/results/imputado/extracted/snps_predictdb_hg19_tofilter_2.txt \
$path/results/imputado/extracted/chr1.dose.vcf.gz \
-o $path/results/imputado/extracted/filtered/chr1.dose_filter.vcf.gz

bcftools view \
-R $path/results/imputado/extracted/hglft_www_1e6f8_5a7fd0_nochr.bed \
$path/results/imputado/extracted/chr1.dose.vcf.gz \
-o $path/results/imputado/extracted/filtered/chr1.dose_filter.vcf.gz

# 3. predicción con posiciones filtradas (nos ahorramos lo de los rsid)
