plink --freq --bfile <input> --out <freq-file>
# seguimos los pasos de https://imputationserver.readthedocs.io/en/latest/prepare-your-data/ 
wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

# ejecutamos plink en wsl para generar el archivo de freq:
"$path/software/plink" --freq --bfile \
"$path/results/plink_data/binary/processed/GSE33528_qc_cr_in" \
--allow-extra-chr \
--out "$path/Repositorio/TFM/results/imputation/freq-file"

# sustituir por nombres de rutas bien más adelante para ejecutar las rutas
cp "$path/results/plink_data/binary/processed/GSE33528_qc_cr_in.bim" /home/miguel/impute_ubuntu
cp "$path/Repositorio/TFM/results/imputation/freq-file"* /home/miguel/impute_ubuntu

#ejecutamos el script de pearl (en el supercomputador - scripts/imputation/cluster_perl_prepare.sh-)para generar el bim actualizado 
# con snps con control de calidad y en formato vcf

# guardamos los vcf en results/imputation/vcf_impute

# luego los comprimimos con bgzip para poder imputar:
for file in $path/results/imputation/vcf_impute/*.vcf; do
    if [ -f "$file" ]; then
        bgzip "$file"
    fi
done

# Finalmente comprobamos ahora el número de SNPs total resultante tras la preparación de datos:
bcftools view -H $path/results/imputation/vcf_impute/old/*.gz | wc -l
bcftools index $path/results/imputation/vcf_impute/old/*.gz

for vcf in $path/results/imputation/vcf_impute/old/*.vcf.gz; do
    bcftools index "$vcf"
done

total_snps=0
for vcf in $path/results/imputation/vcf_impute/old/*.vcf.gz; do
    snp_count=$(bcftools view -H "$vcf" | wc -l)
    total_snps=$((total_snps + snp_count))
done     
# guardamos el número de SNPs resultante para 
echo "$total_snps" > "$path/results/metadata_gsm/SNPs_before_imputation.txt"

# para poderlos imputar tenemos que cambiar el nombre de la columna CHROM del vcf de acorde a hg38 (chr1) 
# para ello hacemos un txt que contenga el mapeo:
for i in {1..22} X Y; do
    echo "$i chr$i" >> chr_conv.txt
done

# y despues hacemos la conversion de todos los archivos:
for i in {1..22} X Y; do
    bcftools annotate --rename-chrs chr_conv.txt GSE33528_qc_cr_in-updated-chr${i}.vcf -Oz -o impute-chr${i}.vcf.gz
done


