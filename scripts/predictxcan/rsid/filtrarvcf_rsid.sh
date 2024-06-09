
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


############# ahora filtraremos los rsid con las variantes de predictdb (que obtuvimos en el script process_db_rsid):

for chr in {1..5} ; do
bcftools view -i ID=@$path/data/predictXcan/elastic_net_models/all_rsids.txt \
$path/results/imputado/extracted/rsid/chr${chr}.dose_rsid.vcf.gz \
-o $path/results/imputado/extracted/filtered/chr${chr}.dose_rsid_filter.vcf.gz
echo "${chr} filtrado!"
done

# uno a uno
bcftools view -i ID=@$path/data/predictXcan/elastic_net_models/all_rsids.txt \
$path/results/imputado/extracted/rsid/chr1.dose_rsid.vcf.gz \
-o $path/results/imputado/extracted/filtered/chr1.dose_rsid_filter2.vcf.gz

############# finalmente concatenamos todos los cromosomas:







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




#####################################

##########nuevo PLANTEAMIENTO 1

#1. liftover la bd de predictdb a hg19, (si no se puede filtrar por rsid)

#2 seleccionar con view solo esas posicioanes en todos los cromosomas

# 3 rsid anotar

# 4 unir todo 

#5predictdb

