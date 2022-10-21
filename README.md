# Inflammation in chronic pain and depression


## Phenotypic directory:

- UKB_CP_MDD_inflam_prep.R: preprocessing of inflammation (CRP and glycA), chronic pain, depression and covariate data in UKB.
- UKB_chronic_pain_depression_CRP_phenotypic_analysis.R: phenotypic association analysis of CRP with chronic pain, depression and comorbid chronic pain and depression in UKB.
- UKB_chronic_pain_depression_glycA_phenotypic_association.R: phenotypic association analysis of glycA with chronic pain, depression and comorbid chronic pain and depression in UKB.
- UKB_chronic_pain_depression_inflammation_phenotypic_plots.R: Plot results of phenoypic association in panels to reduce space.


## MR (Mendelian Randomisation) directory:

- CRP_glycA_CP_MDD_MR.Rmd: Two-sample MR analysis looking at CRP/glycA <-> chronic pain/depression
- olinkCVD_panel_CP_MDD_MR.Rmd: Two-sample MR analysis looking at olink CVD panel <-> chronic pain/depression
- inflamamtion_chronic_pain_depression_MR_resutls_present.Rmd: plot results of two-sample MR analysis

## Colocalisation directory:

- 1_get_LD_SNP_list.R: Get list of SNPs within +/- 500kb of _CD40_ transcript location
- 2_create_LD_matrix.sh: Create LD matrix for SNPs in SNPlist
- 3_harmonise_GWAS_1000g.R: Harmonise GWAS summary statistics
- 4_susie_coloc.R: Carry out fine-mapping and colocalisation analysis
