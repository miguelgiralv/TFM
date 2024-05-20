/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final --exclude /home/mgiralv/tfm/imputation/Exclude-GSE33528_final-1000G.txt --make-bed --out /home/mgiralv/tfm/imputation/TEMP1
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/TEMP1 --update-map /home/mgiralv/tfm/imputation/Chromosome-GSE33528_final-1000G.txt --update-chr --make-bed --out /home/mgiralv/tfm/imputation/TEMP2
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/TEMP2 --update-map /home/mgiralv/tfm/imputation/Position-GSE33528_final-1000G.txt --make-bed --out /home/mgiralv/tfm/imputation/TEMP3
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/TEMP3 --flip /home/mgiralv/tfm/imputation/Strand-Flip-GSE33528_final-1000G.txt --make-bed --out /home/mgiralv/tfm/imputation/TEMP4
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/TEMP4 --a2-allele /home/mgiralv/tfm/imputation/Force-Allele1-GSE33528_final-1000G.txt --make-bed --out /home/mgiralv/tfm/imputation/GSE33528_final-updated
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 1 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr1
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 1 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr1
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 2 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr2
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 2 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr2
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 3 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr3
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 3 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr3
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 4 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr4
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 4 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr4
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 5 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr5
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 5 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr5
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 6 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr6
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 6 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr6
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 7 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr7
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 7 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr7
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 8 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr8
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 8 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr8
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 9 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr9
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 9 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr9
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 10 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr10
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 10 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr10
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 11 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr11
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 11 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr11
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 12 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr12
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 12 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr12
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 13 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr13
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 13 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr13
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 14 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr14
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 14 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr14
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 15 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr15
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 15 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr15
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 16 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr16
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 16 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr16
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 17 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr17
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 17 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr17
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 18 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr18
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 18 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr18
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 19 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr19
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 19 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr19
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 20 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr20
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 20 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr20
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 21 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr21
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 21 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr21
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 22 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr22
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 22 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr22
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --make-bed --chr 23 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr23
/home/mgiralv/tfm/imputation/plink --bfile /home/mgiralv/tfm/imputation/GSE33528_final-updated --real-ref-alleles --recode vcf --chr 23 --out /home/mgiralv/tfm/imputation/GSE33528_final-updated-chr23
rm TEMP*
