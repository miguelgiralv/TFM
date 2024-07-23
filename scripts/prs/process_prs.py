import os
import pandas as pd

path = "C:/Users/Miguel/Documents/UNIVERSIDAD/6_MASTER_BIOINFORMATICA/TFM/Repositorio/TFM"

# procesamos los archivos de score para que en vez de tener los rsid tengan las posiciones en la columna de id
prs_file = pd.read_csv(f"{path}/data/PRS/PGS004146_hmPOS_GRCh37_2.txt", sep="\t")
                         
# para PGS004146, nos quedaremos con chr:pos:ref:alt 
prs_file = pd.read_csv(f"{path}/data/PRS/PGS004146_hmPOS_GRCh37_2.txt", sep="\t")
prs_file['Concatenated'] = prs_file['chr_name'].astype(str) + ':' + \
                     prs_file['chr_position'].astype(str) + ':' + \
                     prs_file['other_allele'] + ':' + \
                     prs_file['effect_allele']               
                     
columns_to_keep = ['Concatenated', 'chr_name', 'chr_position', 'effect_allele', 'other_allele',
       'effect_weight', 'hm_source', 'hm_rsID', 'hm_chr', 'hm_pos',
       'hm_inferOtherAllele']

prs_file2 = prs_file[columns_to_keep]
 
prs_file2.to_csv(f"{path}/data/PRS/PGS004146_hmPOS_GRCh37_3.txt",sep="\t", index=False)

# para PGS000054, nos quedaremos con chr:pos 

prs_file3 = pd.read_csv(f"{path}/data/PRS/PGS000054_hmPOS_GRCh37_2.txt", sep="\t")
                        
prs_file3['Concatenated'] = prs_file3['chr_name'].astype(str) + ':' + \
                        prs_file3['chr_position'].astype(str) 
                        
columns_to_keep = ['Concatenated', 'chr_name', 'chr_position', 'effect_allele', 'effect_weight',
       'OR', 'hm_source', 'hm_rsID', 'hm_chr', 'hm_pos', 'hm_inferOtherAllele',
       'hm_match_chr', 'hm_match_pos']

prs_file4 = prs_file3[columns_to_keep]
                      
prs_file4.to_csv(f"{path}/data/PRS/PGS000054_hmPOS_GRCh37_5.txt",sep="\t", index=False)
