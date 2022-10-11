#install.packages("snpsettest")
library(snpsettest)
library(dplyr)

## Load in exposure - TNFRSF5
exposure_CD40 <- read.table("~/Desktop/PhD/resources/summary_statistics/raw_data/olink_CVD/TNFRSF5.txt", header = TRUE)

## Load in outcome - MDD 
outcome_MDD <-  read.table("~/Desktop/PhD/resources/summary_statistics/processed_data/MDD.txt", header = TRUE)

## Load in outcome - MCP
outcome_MCP <-  read.table("~/Desktop/PhD/resources/summary_statistics/raw_data/ChronicPain.stats", header = TRUE)

## Get SNPs within 500kb of CD40 gene
## Get cis region of exposure protein
## Chromosome 20
bp_range = (46118316 - 500000):(46129863  + 500000)

## Get SNPs in cis region in exposure and outcome GWAS
CD40_snps <- exposure_CD40$hm_rsid[(exposure_CD40$hm_chrom == 20 & exposure_CD40$hm_pos %in% bp_range)]

## Get SNPs in in both exposure and outcome GWAS
MDD_snps <- CD40_snps[CD40_snps %in% outcome_MDD$ID]
MCP_snps <- CD40_snps[CD40_snps %in% outcome_MCP$SNP]

## Change column names for consistency and limit to CD40 snps

exposure_sCD40_QC <- exposure_CD40 %>%
  filter(hm_rsid  %in% CD40_snps) %>%
  dplyr::select(hm_rsid, hm_other_allele, hm_effect_allele, hm_effect_allele_frequency, hm_beta, standard_error, p_value, n) %>%
  dplyr::rename(id = hm_rsid,
                beta = hm_beta,
                freq = hm_effect_allele_frequency,
                pvalue = p_value,
                A1 = hm_effect_allele,
                A2 = hm_other_allele) %>%
  mutate(varbeta = standard_error^2, zscore = beta/standard_error)

## Calculate allele frequency

outcome_MDD_QC <- outcome_MDD %>%
  filter(ID  %in% MDD_snps) %>%
  dplyr::select(ID, A1, A2, freq, BETA, SE, PVAL, CHROM, POS) %>%
  dplyr::rename(id = ID,
                beta = BETA,
                pvalue = PVAL,
                chr = CHROM,
                pos = POS) %>%
  mutate(varbeta = SE^2, zscore = beta/SE)


outcome_MCP_QC <- outcome_MCP %>%
  filter(SNP  %in% MCP_snps) %>%
  dplyr::select(SNP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P, CHR, BP) %>%
  dplyr::rename(id = SNP,
                beta = BETA,
                freq = A1FREQ,
                pvalue =P,
                A1 = ALLELE1,
                A2 = ALLELE0,
                chr = CHR,
                pos = BP) %>%
  mutate(varbeta = SE^2, zscore = beta/SE)

## combine SNP, chr and bp columns of MDD and MCP to apply to SNPs in CD40
outcome_snp_info <- rbind(outcome_MCP_QC[,c("id", "chr", "pos")], outcome_MDD_QC[,c("id", "chr", "pos")])

## Remove duplicates
outcome_snp_info <- outcome_snp_info %>% 
  filter(duplicated(id) == FALSE)

## Add pos and chr to CD40 (hg18)
exposure_sCD40_QC <- inner_join(exposure_sCD40_QC, outcome_snp_info[,c("id", "chr", "pos")], by = "id")

## Remove rows with missing data
outcome_MDD_QC <- na.omit(outcome_MDD_QC)
outcome_MCP_QC <- na.omit(outcome_MCP_QC)
exposure_sCD40_QC <- na.omit(exposure_sCD40_QC)

## Remove duplicated SNPs - keep first instance
outcome_MDD_QC <- outcome_MDD_QC %>% 
  filter(duplicated(id) == FALSE)

outcome_MCP_QC <- outcome_MCP_QC %>% 
  filter(duplicated(id) == FALSE)

exposure_sCD40_QC <- exposure_sCD40_QC %>% 
  filter(duplicated(id) == FALSE)

  
## Load in BED reference file 
x <- read_reference_bed("~/Desktop/PhD/resources/LD_reference/all_phase3.bed", verbose = FALSE)

## Harmonise MDD sumstats
hsumstats_MDD <- harmonize_sumstats(outcome_MDD_QC, x, return_indice = F, check_strand_flip = TRUE)
temp <- left_join(hsumstats_MDD, outcome_MDD_QC, by = "id")
flipped_rsids_MDD <- temp$id[temp$A1.x != temp$A1.y]

## reverse beta and frequency for flipped alleles
outcome_MDD_QC_harmonised <- outcome_MDD_QC[outcome_MDD_QC$id %in% hsumstats_MDD$id,]
outcome_MDD_QC_harmonised[outcome_MDD_QC$id %in% flipped_rsids_MDD, "beta"] <- outcome_MDD_QC[outcome_MDD_QC$id %in% flipped_rsids_MDD, "beta"] * -1
outcome_MDD_QC_harmonised[outcome_MDD_QC$id %in% flipped_rsids_MDD, "zscore"] <- outcome_MDD_QC[outcome_MDD_QC$id %in% flipped_rsids_MDD, "zscore"] * -1
outcome_MDD_QC_harmonised[outcome_MDD_QC$id %in% flipped_rsids_MDD, "freq"] <- 1 - outcome_MDD_QC[outcome_MDD_QC$id %in% flipped_rsids_MDD, "freq"]
outcome_MDD_QC_harmonised[outcome_MDD_QC$id %in% flipped_rsids_MDD, "A1"] <- outcome_MDD_QC[outcome_MDD_QC$id %in% flipped_rsids_MDD, "A2"]
outcome_MDD_QC_harmonised[outcome_MDD_QC$id %in% flipped_rsids_MDD, "A2"] <- outcome_MDD_QC[outcome_MDD_QC$id %in% flipped_rsids_MDD, "A1"]

## Harmonise MCP sumstats
hsumstats_MCP <- harmonize_sumstats(outcome_MCP_QC, x, return_indice = F, check_strand_flip = TRUE)
temp <- inner_join(hsumstats_MCP, outcome_MCP_QC, by = "id")
flipped_rsids_MCP <- temp$id[temp$A1.x != temp$A1.y]

## reverse beta and frequency for flipped alleles
outcome_MCP_QC_harmonised <- outcome_MCP_QC[outcome_MCP_QC$id %in% hsumstats_MCP$id,]
outcome_MCP_QC_harmonised[, "beta"] <- outcome_MCP_QC[outcome_MCP_QC$id %in% flipped_rsids_MCP, "beta"] * -1
outcome_MCP_QC_harmonised[, "zscore"] <- outcome_MCP_QC[outcome_MCP_QC$id %in% flipped_rsids_MCP, "zscore"] * -1
outcome_MCP_QC_harmonised[, "freq"] <- 1 - outcome_MCP_QC[outcome_MCP_QC$id %in% flipped_rsids_MCP, "freq"]
outcome_MCP_QC_harmonised[, "A1"] <- outcome_MCP_QC[outcome_MCP_QC$id %in% flipped_rsids_MCP, "A2"]
outcome_MCP_QC_harmonised[, "A2"] <- outcome_MCP_QC[outcome_MCP_QC$id %in% flipped_rsids_MCP, "A1"]

## Harmonise sCD40 sumstats
hsumstats_sCD40 <- harmonize_sumstats(exposure_sCD40_QC, x, return_indice = F, check_strand_flip = TRUE)
temp <- left_join(hsumstats_sCD40, exposure_sCD40_QC, by = "id")
flipped_rsids_sCD40 <- temp$id[temp$A1.x != temp$A1.y]

## No SNPs needing harmonization!


## Save harmonised sumstats
write.table(outcome_MDD_QC_harmonised, '~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/MDD_harmonised.tsv', sep = '\t', row.names = F, quote = F)
write.table(outcome_MCP_QC_harmonised, '~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/MCP_harmonised.tsv', sep = '\t', row.names = F, quote = F)
write.table(exposure_sCD40_QC, '~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/sCD40_harmonised.tsv', sep = '\t', row.names = F, quote = F)
