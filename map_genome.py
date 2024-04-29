
#################################################################################
# mapeo de posiciones desde hg18 a hg38

# mapeamos las coordenadas de hg18 a hg38
# descargamos la cadena de https://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/?C=M;O=A
# Probamos con el liftover de GATK: https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard  
#primero generamos el diccionario

import subprocess

liftover_dictionary= ["java", "-jar", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/picard.jar", "CreateSequenceDictionary",
"REFERENCE=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/GCF_000001405.26_GRCh38_genomic.fna",
"OUTPUT=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/GCF_000001405.26_GRCh38_genomic.dict"]
subprocess.run(liftover_dictionary)

#despues ejecutamos liftovervcf
liftover_command = ["java","-Xmx8g","-jar" ,"C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/picard.jar",
"LiftoverVcf",
"I=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_test.vcf",
"O=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/lifted_over.vcf",
"CHAIN=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/hg18ToHg38.over.chain.gz",
"REJECT=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/rejected_variants.vcf",
"R=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/GCF_000001405.26_GRCh38_genomic.fna"] 

#ejecutamos map_to_UCSC_bed.py:

map_to_bed=["python","C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/map_to_UCSC_bed_python3.py", "-m", 
"C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528.map",
"-o", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/GSE33528_UCSC"]

subprocess.run(map_to_bed)

# Lo ejecutamos para obtener los nuevos archivos binarios excluyendo a las muestras que no pasaron la calidad
subprocess.run(liftover_command)

#java -jar -Xmx8g picard.jar LiftoverVcf I=GSE33528_test.vcf O=lifted_over.vcf CHAIN=hg18ToHg38.over.chain.gz REJECT=rejected_variants.vcf R=GCF_000001405.26_GRCh38_genomic.fna

subprocess.run("java -Xmx8g -jar C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/picard.jar \
C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/software/LiftoverVcf \
I=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/GSE33528_test.vcf \
O=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/lifted_over.vcf \
CHAIN=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/hg18ToHg38.over.chain.gz \
REJECT=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/binary/processed/rejected_variants.vcf \
R=C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/coordinate_map/GCF_000001405.26_GRCh38_genomic.fna", shell=True)

############################################################

# intento mapear de .map a bed usc:

# https://genome.sph.umich.edu/wiki/LiftOver#Method_1
# copio parte del script de https://genome.sph.umich.edu/wiki/LiftMap.py

def myopen(fn):
    import gzip
    try:
        with gzip.open(fn, 'rb') as h:
            ln = h.read(2)  # read arbitrary bytes so check if @param fn is a gzipped file
    except:
        # cannot read in gzip format
        return open(fn, 'r')
    return gzip.open(fn, 'rt')


def map2bed(fin, fout):
    with open(fout, 'w') as fo:
        for ln in myopen(fin):
            chrom, rs, mdist, pos = ln.split()
            chrom = 'chr' + chrom
            pos = int(pos)
            fo.write('%s\t%d\t%d\t%s\n' % (chrom, pos-1, pos, rs))
    return True
## generamos el archivo bed
map2bed("C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528.map", "C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/results/plink_data/classic/GSE33528_out_ucsc.bed")
 
# ejecutamos liftover en wsl:
#  ./liftOver GSE33528_out_ucsc.bed hg18ToHg38.over.chain.gz output.bed unlifted.bed