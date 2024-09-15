path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

# 1.filtrado con las posiciones de predictdb
# filtramos todos los cromosomas por partes en varias consolas de wsl para ser más rápidos: 
for chr in {1..5} ; do
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

# 2.Concatenamos todos los cromosomas para tener un único vcf:
# hay que indexar antes:
for chr in {1..22} X; do
    bcftools index $path/results/imputado/extracted/filtered/mashr/chr${chr}.dose_filter_mashr.vcf.gz
    echo "${chr} indexado!"
done
#concatenamos
bcftools concat -Oz -o $path/results/imputado/extracted/filtered/mashr/fullgenome_mashr.vcf.gz \
    $path/results/imputado/extracted/filtered/mashr/chr{1..22}.dose_filter_mashr.vcf.gz \
    $path/results/imputado/extracted/filtered/mashr/chrX.dose_filter_mashr.vcf.gz \

#contamos el numero de SNPs tras filtrar:
bcftools view -v snps $path/results/imputado/extracted/filtered/mashr/fullgenome_mashr2.vcf.gz | grep -v "^#" | wc -l
# 216665 SNPs