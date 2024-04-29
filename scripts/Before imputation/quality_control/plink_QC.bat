REM Ejecutamos plink para eliminar las muestras que fallaron el test de calidad y que hemos almacenado en qc_failed_ids.txt

"C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink.exe" ^
"--bfile" "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/raw/output_GSE33528" ^
"--remove" "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/metadata_gsm/qc_failed_ids.txt" ^ 
"--make-bed" "--out" "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc"

REM Seleccionamos también los individuos con buena calidad en base al call-rate de variantes (call rate de 0.95) (comando --geno)

"C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink.exe" ^
"--bfile" "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc" ^
"--geno" "0.05" ^
"--out" "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc_cr" ^
"--make-bed"

REM Seleccionamos también los individuos con buena calidad en base al call-rate de individuos (comando --mind)

"C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/plink.exe" ^
"--bfile" "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc_cr" ^
"--mind" "0.05" ^
"--out" "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_qc_cr_in" ^
"--make-bed"




