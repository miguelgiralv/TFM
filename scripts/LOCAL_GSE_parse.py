import pandas as pd
import GEOparse
import pandas_plink
import xarray
import sys

# codigo original modificado:
# https://gist.github.com/CholoTook/60968e3ab6d90cb8fd19be55a25592f1

def get_gpl(gpl_accession='GPL14932'):
    gpl = GEOparse.get_GEO(gpl_accession,
                           destdir="C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/GEO")
    # Para GPL14932 el identificador de SNP es SNP_ID, no SNP
    gpl.table['a0'] = [x[1] for x in gpl.table.SNP_ID]
    gpl.table['a1'] = [x[3] for x in gpl.table.SNP_ID]
    return gpl


def geo_gsm_to_plink(gse_accession, gpl):
    gse = GEOparse.get_GEO(gse_accession,
                           destdir='C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/GEO/samples')   
    for gsm_name, gsm in gse.gsms.items():    
        # Perform the merge operation
        m = gsm.table.merge(gpl.table, left_on='ID_REF', right_on='ID')
        # Carefully build a DataArray object to pass into pandas-plink for writing
        genotype_data = xarray.DataArray(
            data = [m["Theta"]],
            dims = ["sample", "variant"],
            coords = dict(
                iid   = ("sample", [gsm_name]),
                snp   = ("variant", m.ID),
                chrom = ("variant", m.Chr),
                pos   = ("variant", m.MapInfo),
                a0    = ("variant", m.a0),
                a1    = ("variant", m.a1)
            )
        )
        # Write out Plink format files
        pandas_plink.write_plink1_bin(
            genotype_data, f"C:/Users/Miguel/Documents/UNIVERSIDAD/6 MASTER BIOINFORMATICA/TFM/Repositorio/TFM/data/plink_format/plink{gsm_name}.bed"
        )
        
        # On the command line, use plink --recode to get PED or VCF files!

if __name__ == '__main__':
    geo_gsm_to_plink("GSE33528", get_gpl())


