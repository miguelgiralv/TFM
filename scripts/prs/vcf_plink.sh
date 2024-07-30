
path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

# los cromosomas imputados los convertimos a plink
for chr in {1..22} X; do
    $path/software/plink --vcf $path/results/imputado/extracted/chr${chr}.dose.vcf.gz --make-bed --out $path/results/imputado/extracted/plink/imputed_data_chr${chr}
done

######### a partir de aqui lo hice en el supercomputador (probar en local a ver si funciona)

# convertimos nuestros archivos de plink a plink2 para poder unir todo el genoma (tiene variantes muy largas que plink 1.9 no soporta):
for chr in {1..22} X; do
$path/software/plink2 --bfile $path/results/imputado/extracted/plink/imputed_data_chr${chr} --make-pgen --out $path/results/imputado/extracted/plink/pgen/imputed_data_chr${chr}
echo "done chr${chr}"
done
$path/software/plink2 --bfile $path/results/imputado/extracted/plink/imputed_data_chr1 --make-pgen --out $path/results/imputado/extracted/plink/pgen/imputed_data_chr1

# creamos el archivo chr_list.txt con la ruta de todos los pfiles y luego los unimos con plink2
for chr in {1..22} X; do
echo $path/results/imputado/extracted/plink/pgen/imputed_data_chr${chr} >>$path/results/imputado/extracted/plink/pgen/chr_list.txt
done
$path/software/plink2 --pmerge-list $path/results/imputado/extracted/plink/pgen/chr_list.txt --out $path/results/PRS/full_genome

# volvemos a convertirlo a plink1 (y así aprovechamos el archivo fam que teníamos antes con informacion fenotipica en el plink original) 
$path/software/plink2 --pfile $path/results/PRS/full_genome --make-bed --out $path/results/PRS/fullgenome_plink


# calculamos los prs de PGS000054 (hemos sustituido los rsid por pos:chr:ref:alt en un script de python "process_prs"):
$path/software/plink2 --bfile $path/results/PRS/fullgenome_plink --score $path/data/PRS/PGS000054_hmPOS_GRCh37_5.txt 1 4 5 header --out  $path/results/PRS/PGS000054_hmPOS_GRCh37/output_PGS000054

# calculamos los prs de PGS004146 con el chr posiciones y alelos (tambien transformado con el script de python)
$path/software/plink2 --bfile $path/results/PRS/fullgenome_plink --score $path/data/PRS/PGS004146_hmPOS_GRCh37_4.txt 1 4 6 header --out  $path/results/PRS/PGS004146_hmPOS_GRCh37/output_PGS004146


