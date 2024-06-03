
# primero debemos anotar los vcf con sus rsid
path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

bgzip -c $path/results/imputation/HRC.r1-1.GRCh37.wgs.mac5_2.sites.tab > $path/results/imputation/HRC.r1-1.GRCh37.wgs.mac5_2.sites.tab.gz

tabix -s 1 -b 2 -e 2 $path/results/imputation/HRC.r1-1.GRCh37.wgs.mac5_2.sites.tab.gz


bcftools annotate -c CHROM,POS,~ID,REF,ALT -a $path/results/imputation/HRC.r1-1.GRCh37.wgs.mac5_2.sites.tab.gz \
-o $path/results/imputado/extracted/rsid/chr22_rsid.vcf.gz $path/results/imputado/extracted/chr22.dose.vcf.gz

bcftools annotate -c CHROM,POS,ID,REF,ALT -a $path/results/imputation/HRC.r1-1.GRCh37.wgs.mac5_2.sites.tab.gz \
-o $path/results/imputado/extracted/rsid/chr22_rsid.vcf.gz $path/results/imputado/extracted/chr22.dose.vcf.gz


bcftools annotate -c CHROM,POS,ID -a $path/results/imputation/HRC.r1-1.GRCh37.wgs.mac5_2.sites.tab.gz \
-o $path/results/imputado/extracted/rsid/chr1_test_rsid.vcf $path/results/imputado/extracted/rsid/chr1_test.vcf.gz

bcftools annotate -c ID -a $path/results/imputation/hg19_vcf/GCF_000001405.25.gz \
-o $path/results/imputado/extracted/rsid/chr1_test_rsid.vcf $path/results/imputado/extracted/rsid/chr1_test.vcf.gz

https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/#table-1


bcftools view -v indels $path/results/imputation/hg19_vcf/GCF_000001405.25.gz

#anotado hasta el cromosoma 4

#filtrar la lista de rsids:
bcftools view -R $path/data/predictXcan/elastic_net_models/all_rsids.txt $path/results/imputado/extracted/rsid/output.vcf.gz -o $path/results/imputado/extracted/filtered/chr22.dose_preditdb.vcf.gz

bcftools view -R $path/data/predictXcan/elastic_net_models/all_rsids.txt $path/results/imputado/extracted/rsid/output.vcf.gz 
-o $path/results/imputado/extracted/filtered/chr22.dose_preditdb.vcf.gz

bcftools view -r 1:993810 $path/results/imputado/extracted/chr1.dose.vcf.gz -o $path/results/imputado/extracted/rsid/chr1_test.vcf.gz



# hacemos un bucle for para anotar todos los vcf con sus rsid
for chr in {1..22} X; do
    bcftools annotate -c CHROM,POS,~ID,REF,ALT -a $path/results/imputation/HRC.r1-1.GRCh37.wgs.mac5_2.sites.tab.gz -o $path/results/imputado/extracted/rsid/chr${chr}_rsid.dose.vcf.gz $path/results/imputado/extracted/chr${chr}.dose.vcf.gz
    echo "${chr} anotado!"
done
# 1 y 2 anotado

# tras ello podemos filtrar los vcf de cada cromosoma guardando solo los rsid de predictdb

bcftools view -R $path/data/predictXcan/elastic_net_models/all_rsids.txt $path/results/imputado/extracted/rsid/chr1_rsid.dose.vcf.gz -o $path/results/imputado/extracted/filtered/chr1.dose_preditdb.vcf.gz

#finalmente, juntamos los vcf procesados en un unico vcf que pasaramos por predictscan