
path="/mnt/c/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

# calculamos los prs de PGS000054 (hemos sustituido los rsid por pos:chr:ref:alt en un script de python "process_prs"):
$path/software/plink2 --bfile $path/results/PRS/fullgenome_plink --score $path/data/PRS/PGS000054_hmPOS_GRCh37_5.txt 1 4 5 header --out  $path/results/PRS/PGS000054_hmPOS_GRCh37/output_PGS000054

# calculamos los prs de PGS004146 con el chr posiciones y alelos (tambien transformado con el script de python)
$path/software/plink2 --bfile $path/results/PRS/fullgenome_plink --score $path/data/PRS/PGS004146_hmPOS_GRCh37_4.txt 1 4 6 header --out  $path/results/PRS/PGS004146_hmPOS_GRCh37/output_PGS004146


