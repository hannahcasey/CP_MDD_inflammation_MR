
/Users/hannahcasey/Downloads/plink2 \
      --bgen /Volumes/GenScotDepression/data/ukb/genetics/impv3/ukb_imp_chr20_v3.bgen 'ref-first' \
      --sample /Volumes/GenScotDepression/data/ukb/genetics/impv3/ukb4844_imp_chr20_v3_s487395.sample \
      --extract ~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/CD40_common_snps.txt \
      --make-bed \
      --out ~/Desktop/PhD/resources/LD_reference/UKB_CD40


/Users/hannahcasey/Downloads/plink_mac_20220402/plink \
      --bfile ~/Desktop/PhD/resources/LD_reference/UKB_CD40\
      --r square spaces \
      --keep-allele-order \
      --write-snplist \
      --out ~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/LD_matrix_spaces



/Users/hannahcasey/Downloads/plink_mac_20220402/plink \
      --bfile ~/Desktop/PhD/resources/LD_reference/all_phase3 \
      --keep ~/Desktop/PhD/resources/LD_reference/all_phase3_eur.fam \
      --allow-extra-chr \
      --keep-allele-order \
      --extract ~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/CD40_common_snps.txt \
      --make-bed \
      --out ~/Desktop/PhD/resources/LD_reference/eur_CD40


/Users/hannahcasey/Downloads/plink_mac_20220402/plink \
      --bfile ~/Desktop/PhD/resources/LD_reference/eur_CD40\
      --r square spaces \
      --keep-allele-order \
      --write-snplist \
      --out ~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/LD_matrix

