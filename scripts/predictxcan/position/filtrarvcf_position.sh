
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
    -o $path/results/imputado/extracted/rsid/chr1.dose_rsid.vcf.gz \
    $path/results/imputado/extracted/chr1.dose.vcf.gz

# vemos los resultados (funciona)
bcftools view $path/results/imputation/1000genomas/chr1.dose_rsid.vcf.gz | less

# ahora hacemos un bucle for para anotar todos los cromosomas. Como tarda mucho haremos 4 bucles
for chr in {1..5} ; do
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



############# ahora filtramos los rsid con las variantes de predictdb:

for chr in {1..5} ; do
    bcftools view \
    -R $path/data/predictXcan/elastic_net_models/all_rsids.txt \
    -o $path/results/imputado/extracted/rsid/chr${chr}.dose_filter_rsid.vcf.gz \
    $path/results/imputado/extracted/rsid/chr${chr}.dose_rsid.vcf.gz
    echo "${chr} filtrado!"
done



# visualizamos el nuevo 
bcftools view $path/results/imputado/extracted/rsid/chr1.dose_test_new_1000G.vcf.gz | less



############################# No revisado

#finalmente, juntamos los vcf procesados en un unico vcf que pasaramos por predictscan