path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

# Le quitamos el chr a all_positions_mashr_hg19 poder filtrar con bcftools:
input_file = f"{path}/results/imputado/extracted/all_positions_mashr_hg19.txt" 
output_file = f"{path}/results/imputado/extracted/all_positions_mashr_hg19_nochr.txt" 

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        parts = line.strip().split('\t')  
        chromosome = parts[0].replace('chr', '') 
        start = parts[1]
        end = parts[2]
        outfile.write(f"{chromosome}\t{start}\t{end}\n") 

# y luego filtraremos los vcf con el  (filtrarvcf_position.sh)



 
