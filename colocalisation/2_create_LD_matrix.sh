## Limit genotype data to SNPs located in CD40 region

/Users/hannahcasey/Downloads/plink_mac_20220402/plink \
      --bfile ~/Desktop/PhD/resources/LD_reference/UKB_CD40 \
      --keep ~/Desktop/PhD/resources/LD_reference/UKB_CD40.fam \
      --allow-extra-chr \
      --keep-allele-order \
      --extract ~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/CD40_common_snps.txt \
      --make-bed \
      --out ~/Desktop/PhD/resources/LD_reference/UKB_CD40

## Get LD matrix
/Users/hannahcasey/Downloads/plink_mac_20220402/plink \
      --bfile ~/Desktop/PhD/resources/LD_reference/UKB_CD40\
      --r square spaces \
      --keep-allele-order \
      --write-snplist \
      --out ~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/UKB_CD40_LD

