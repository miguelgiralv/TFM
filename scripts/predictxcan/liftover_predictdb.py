path="C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"



# procesamos el archivo con los snps de predictdeb que generamos antes (donde?) (generar para TODOS los tejidos)
# para poder hacer el liftover a hg19

input_file =  f"{path}/results/imputado/extracted/snps_predictdb.txt"
bed_file = f"{path}/results/imputado/extracted/snps_predict_chr.bed"
with open(input_file, 'r') as infile, open(bed_file, 'w') as outfile:
    for line in infile:
        chrom, pos = line.strip().split(':')
        start = int(pos) - 1
        end = pos
        outfile.write(f"chr{chrom}\t{start}\t{end}\n")


# una vez hemos hecho el liftover, tnedremos que procesarlo de nuevo para poder filtrar 
#con bcftools con view (formato 10:4459595)

input_file =  f"{path}/results/imputado/extracted/hglft_www_1e6f8_5a7fd0.bed"
output_file =  f"{path}/results/imputado/extracted/hglft_www_1e6f8_5a7fd0_nochr.bed"

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        parts = line.strip().split('\t')  # Split the line by tabs
        chromosome = parts[0].replace('chr', '')  # Remove the "chr" part
        start, end = parts[1], parts[2]
        outfile.write(f"{chromosome}\t{start}\t{end}\n")  # Write the modified line to the output file



input_file = f"{path}/results/imputado/extracted/snps_predictdb_hg19_tofilter.txt"
output_file =  f"{path}/results/imputado/extracted/snps_predictdb_hg19_tofilter_2.txt"

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        chrom_pos = line.strip().split(':')  # Split each line by colon (:)
        modified_line = '\t'.join(chrom_pos)  # Join the elements using tab as the delimiter
        outfile.write(modified_line + '\n')  # Write the modified line to the output file

# con esto ya podemos filtrar los cromosomas y luego predecir con predictxcan

