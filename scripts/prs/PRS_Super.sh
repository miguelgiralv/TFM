

salloc -J prs_join -p day -N 2 -c 4


# convertimos nuestros archivos de plink a plink2 para poder unir todo el genoma (tiene variantes muy largas que plink 1.9 no soporta):
for chr in {1..22} X; do
./plink2 --bfile "imputed_data_chr${chr}" --make-pgen --out "pgen/imputed_data_chr${chr}"
echo "done chr${chr}"
done
# creamos el archivo chr_list.txt con la ruta de todos los pfiles y luego los unimos con plink2
./plink2 --pmerge-list chr_list.txt --out full_genome
# volvemos a convertirlo a plink1 (aprovechamos el archivo fam que ya teníamos antes)
./plink2 --pfile full_genome --make-bed --out fullgenome_plink

# para calcular los prs de PGS000054 no tenemos los dos alelos en el score, así que procesamos el id para que solo sea el chr:pos
./plink2 --bfile fullgenome_plink --set-all-var-ids @:# --make-bed --out test/fullgenome_plink_2

# calculamos los prs de PGS004146 con el chr posiciones y alelos
./plink2 --bfile fullgenome_plink --score PGS004146_hmPOS_GRCh37_3.txt 1 4 6 header --out PGS004146_hmPOS_GRCh37/output_PGS004146_2

# calculamos los prs de PGS000054 con el chr y posiciones sin los alelos
./plink2 --bfile test/fullgenome_plink_2 --score PGS000054_hmPOS_GRCh37_5.txt 1 4 5 header --out PGS000054_hmPOS_GRCh37/output_PGS000054


