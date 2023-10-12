library(dplyr)
library(data.table)

## Prep GWAS summary statistics
MDD_GWAS_raw <- read.table("~/Desktop/PhD/resources/summary_statistics/raw_data/pgc-mdd2022-full-eur-v3.49.24.11.pgc", header = T) #GRCh37/hg19
MCP_GWAS_raw <- read.table("~/Desktop/PhD/resources/summary_statistics/raw_data/ChronicPain.stats", header = T) #GRCh37/hg19
sCD40_GWAS_raw <- read.table("~/Desktop/PhD/resources/summary_statistics/raw_data/sCD40.txt", header = T) #GRCh38/hg38

## Get SNPs within 500kb of CD40 gene
## 
bp_range_hg38 = (46118316 - 500000):(46129863 + 500000)
bp_range_hg19 = (44746910 - 500000):(44758502 + 500000)


## Get SNPs in cis region in exposure and outcome GWAS
MDD_GWAS_CD40_SNPs <- MDD_GWAS_raw$ID[(MDD_GWAS_raw$CHROM == 20 & MDD_GWAS_raw$POS %in% bp_range_hg19)]
MCP_GWAS_CD40_SNPs <- MCP_GWAS_raw$SNP[(MCP_GWAS_raw$CHR == 20 & MCP_GWAS_raw$BP %in% bp_range_hg19)]
sCD40_GWAS_CD40_SNPs <- sCD40_GWAS_raw$hm_rsid[(sCD40_GWAS_raw$hm_chrom == 20 & sCD40_GWAS_raw$hm_pos %in% bp_range_hg38)]

## Limit GWAS sumstats to CD40 region and columns required for coloc: beta, varbeta, N, sdY, type, MAF, snp

## Change column names for consistency 
sCD40_QC <- sCD40_GWAS_raw %>%
  filter(hm_rsid  %in% sCD40_GWAS_CD40_SNPs) %>%
  dplyr::select(hm_rsid, hm_other_allele, hm_effect_allele, hm_effect_allele_frequency, hm_beta, standard_error, p_value, n, hm_pos, hm_chrom) %>%
  dplyr::rename(id = hm_rsid,
                beta = hm_beta,
                freq = hm_effect_allele_frequency,
                pvalue = p_value,
                A1 = hm_effect_allele,
                A2 = hm_other_allele) %>%
  mutate(varbeta = standard_error^2, zscore = beta/standard_error)

rm(sCD40_GWAS_raw)

MDD_QC <- MDD_GWAS_raw %>%
  filter(ID  %in% MDD_GWAS_CD40_SNPs) %>%
  dplyr::rename(id = ID,
                beta = BETA,
                pvalue = PVAL,
                pos = POS) %>%
  mutate(varbeta = SE^2, zscore = beta/SE, n = (NCAS + NCON), s = (NCAS/(NCAS + NCON)), freq = (((FCAS * NCAS) + (FCON * NCON)) / (NCAS + NCON))) %>%
  dplyr::select(id, A1, A2, freq, beta, SE, pvalue, pos, s, n, zscore, varbeta)

rm(MDD_GWAS_raw)

MCP_QC <- MCP_GWAS_raw %>%
  filter(SNP  %in% MCP_GWAS_CD40_SNPs) %>%
  dplyr::select(SNP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P, BP) %>%
  dplyr::rename(id = SNP,
                beta = BETA,
                pvalue =P,
                A1 = ALLELE1,
                A2 = ALLELE0,
                freq = A1FREQ,
                pos = BP) %>%
  mutate(varbeta = SE^2, zscore = beta/SE)

rm(MCP_GWAS_raw)

## Convert coordinates of sCD40 sumstats from hg38 to hg19
sCD40_hg38 <- paste0("chr",sCD40_QC$hm_chrom, ":", sCD40_QC$hm_pos, "-", sCD40_QC$hm_pos)
write.table(sCD40_hg38, "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/sCD40_hg38.txt", quote = F, row.names = F, col.names = F)
## Convert online: https://genome.ucsc.edu/cgi-bin/hgLiftOver
## Load in converted positions
sCD40_hg19 <- read.table("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/hglft_genome_387be_db1af0.bed")
sCD40_hg19 <- gsub(".*:","", sCD40_hg19$V1)
sCD40_hg19 <- gsub("-.*","", sCD40_hg19)
sCD40_QC$pos <- as.numeric(sCD40_hg19)

## Add chromosome
MDD_QC$chr <- "20"
MCP_QC$chr <- "20"
sCD40_QC$chr <- "20"

## Remove duplicated SNPs - keep first instance
MDD_QC <- MDD_QC %>% 
  filter(duplicated(id) == FALSE)

MCP_QC <- MCP_QC %>% 
  filter(duplicated(id) == FALSE)

sCD40_QC <- sCD40_QC %>% 
  filter(duplicated(id) == FALSE)

## Save sumstats
write.table(MDD_QC, '~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/MDD_CD40.tsv', sep = '\t', row.names = F, quote = F)
write.table(MCP_QC, '~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/MCP_CD40.tsv', sep = '\t', row.names = F, quote = F)
write.table(sCD40_QC, '~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/sCD40_CD40.tsv', sep = '\t', row.names = F, quote = F)
