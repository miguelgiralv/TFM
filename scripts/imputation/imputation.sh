plink --freq --bfile <input> --out <freq-file>
# seguimos los pasos de https://imputationserver.readthedocs.io/en/latest/prepare-your-data/ 
wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.7.zip
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

# ejecutamos plink en wsl para generar el archivo de freq:
"/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/software/plink" --freq --bfile \
"/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_final" \
--allow-extra-chr \
--out "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/imputation/freq-file"

# sustituir por nombres de rutas bien más adelante para ejecutar las rutas
cp "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_final.bim" /home/miguel/impute_ubuntu
cp "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/imputation/freq-file"* /home/miguel/impute_ubuntu

#ejecutamos el script de pearl (en el supercomputador - -)para generar el bim actualizado con snps con control de calidad y en formato vcf:
perl "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/software/HRC-1000G-check-bim-NoReadKey.pl" \
-b "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_final.bim" \
-f "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/imputation/freq-file" \
-r "/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/data/imputar/1000GP_Phase3_combined.legend" -g

# guardamos los vcf en results\imputation\chr, los copiamos a vcf_impute:
cp /mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/imputation/chr/*.vcf /mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/imputation/vcf_impute


# luego los comprimimos con bgzip para poder imputar:
for file in /mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM/results/imputation/vcf_impute/*.vcf; do
    if [ -f "$file" ]; then
        bgzip "$file"
    fi
done

################
#DUDAS?
# para poderlos imputar tenemos que cambiar el nombre de la columna CHROM del vcf de acorde a hg38 (chr1) 
# para ello hacemos un txt que contenga el mapeo:
for i in {1..22} X Y; do
    echo "$i chr$i" >> chr_conv.txt
done

# y despues hacemos la conversion de todos los archivos:
for i in {1..22} X Y; do
    bcftools annotate --rename-chrs chr_conv.txt GSE33528_final-updated-chr${i}.vcf -Oz -o impute-chr${i}.vcf.gz
done


