#install.packages("snpsettest")
library(snpsettest)
library(dplyr)
library(data.table)

## Load in exposure - TNFRSF5 - (GRCh38)
# exposure_CD40 <- read.table("~/Desktop/PhD/resources/summary_statistics/raw_data/olink_CVD/sCD40_harmonised.tsv", header = TRUE)
exposure_CD40 <- read.table("/Users/hannahcasey/Desktop/PhD/resources/summary_statistics/raw_data/sCD40.txt", header = TRUE)

## Load in outcome - MDD (GRCh37)
outcome_MDD <-  read.table("~/Desktop/PhD/resources/summary_statistics/processed_data/MDD.txt", header = TRUE)

## Load in outcome - MCP - (GRCh37)
outcome_MCP <-  read.table("~/Desktop/PhD/resources/summary_statistics/raw_data/ChronicPain.stats", header = TRUE)


## Get SNPs within 500kb of CD40 gene
## Get cis region of exposure protein
## Chromosome 20
bp_range_hg38 = (46118316 - 500000):(46129863 + 500000)
bp_range_hg19 = ( 44746910 - 500000):(44758502 + 500000)


## Get SNPs in cis region in exposure and outcome GWAS
CD40_snps <- exposure_CD40$hm_rsid[(exposure_CD40$hm_chrom == 20 & exposure_CD40$hm_pos %in% bp_range_hg38)]
MDD_snps <- outcome_MDD$ID[(outcome_MDD$CHROM == 20 & outcome_MDD$POS %in% bp_range_hg19)]
MCP_snps <- outcome_MCP$SNP[(outcome_MCP$CHR == 20 & outcome_MCP $BP %in% bp_range_hg19)]

## Change column names for consistency 
exposure_sCD40_QC <- exposure_CD40 %>%
  filter(hm_rsid  %in% CD40_snps) %>%
  dplyr::select(hm_rsid, hm_other_allele, hm_effect_allele, hm_effect_allele_frequency, hm_beta, standard_error, p_value, n, hm_pos, hm_chrom) %>%
  dplyr::rename(id = hm_rsid,
                beta = hm_beta,
                freq = hm_effect_allele_frequency,
                pvalue = p_value,
                A1 = hm_effect_allele,
                A2 = hm_other_allele) %>%
  mutate(varbeta = standard_error^2, zscore = beta/standard_error)


outcome_MDD_QC <- outcome_MDD %>%
  filter(ID  %in% MDD_snps) %>%
  dplyr::rename(id = ID,
                beta = BETA,
                pvalue = PVAL,
                pos = POS) %>%
  mutate(varbeta = SE^2, zscore = beta/SE, n = (NCAS + NCON), s = (NCAS/(NCAS + NCON))) %>%
 dplyr::select(id, A1, A2, freq, beta, SE, pvalue, pos, s, n, zscore, varbeta)
  

outcome_MCP_QC <- outcome_MCP %>%
  filter(SNP  %in% MCP_snps) %>%
  dplyr::select(SNP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P, BP) %>%
  dplyr::rename(id = SNP,
                beta = BETA,
                pvalue =P,
                A1 = ALLELE1,
                A2 = ALLELE0,
                freq = A1FREQ,
                pos = BP) %>%
  mutate(varbeta = SE^2, zscore = beta/SE)

## Convert coordinates of CD40 sumstats from hg38 to hg19
exposure_sCD40_hg38 <- paste0("chr",exposure_sCD40_QC$hm_chrom, ":", exposure_sCD40_QC$hm_pos, "-", exposure_sCD40_QC$hm_pos)
write.table(exposure_sCD40_hg38, "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/exposure_sCD40_hg38.txt", quote = F, row.names = F, col.names = F)
## Convert online: https://genome.ucsc.edu/cgi-bin/hgLiftOver
## Load in converted positions
exposure_sCD40_hg19 <- read.table("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/hglft_genome_387be_db1af0.bed")
exposure_sCD40_hg19 <- gsub(".*:","", exposure_sCD40_hg19$V1)
exposure_sCD40_hg19 <- gsub("-.*","", exposure_sCD40_hg19)
exposure_sCD40_QC$pos <- as.numeric(exposure_sCD40_hg19)

## Add chromosome
outcome_MDD_QC$chr <- "20"
outcome_MCP_QC$chr <- "20"
exposure_sCD40_QC$chr <- "20"


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

## Load in BIM file
bim_file <- as.data.frame(fread("~/Desktop/PhD/resources/LD_reference/UKB_CD40.bim"))
## Add column names
names(bim_file) <- c("chr", "id", "pos_measure", "pos", "A1", "A2")
## Remove SNPs that appear twice - not sure which LD correlation applies
bim_file <- bim_file[!bim_file$id %in% bim_file$id[duplicated(bim_file$id)],]

## Make column combining A1 and A2 in alphabetical order
bim_file$A1A2 <- apply(bim_file[,c("A1", "A2")], 1, function(x) paste(sort(x), collapse = ""))
outcome_MDD_QC$A1A2 <- apply(outcome_MDD_QC[,c("A1", "A2")], 1, function(x) paste(sort(x), collapse = ""))
outcome_MCP_QC$A1A2 <- apply(outcome_MCP_QC[,c("A1", "A2")], 1, function(x) paste(sort(x), collapse = ""))
exposure_sCD40_QC$A1A2 <- apply(exposure_sCD40_QC[,c("A1", "A2")], 1, function(x) paste(sort(x), collapse = ""))

## Merge BIM and sumstats
outcome_MDD_QC_bim <- inner_join(outcome_MDD_QC, bim_file, by = c("id", "pos"))
outcome_MCP_QC_bim <- inner_join(outcome_MCP_QC, bim_file, by = c("id", "pos"))
exposure_CD40_QC_bim <- inner_join(exposure_sCD40_QC, bim_file, by = c("id", "pos"))

## Get SNP ids to flip 
flipped_rsids_MDD <- outcome_MDD_QC_bim$id[outcome_MDD_QC_bim$A1.x != outcome_MDD_QC_bim$A1.y]
flipped_rsids_MCP <- outcome_MCP_QC_bim$id[outcome_MCP_QC_bim$A1.x != outcome_MCP_QC_bim$A1.y]
flipped_rsids_CD40 <- exposure_CD40_QC_bim$id[exposure_CD40_QC_bim$A1.x != exposure_CD40_QC_bim$A1.y]

## Get SNP ids to remove - where flipping will not allign alleles or SNPs are ambiguous
remove_rsids_MDD <- outcome_MDD_QC_bim$id[outcome_MDD_QC_bim$A1A2.x == "AT" | outcome_MDD_QC_bim$A1A2.x == "CG"]
remove_rsids_MDD <- c(remove_rsids_MDD, outcome_MDD_QC_bim$id[outcome_MDD_QC_bim$A1A2.x != outcome_MDD_QC_bim$A1A2.y])

remove_rsids_MCP <- outcome_MCP_QC_bim$id[outcome_MCP_QC_bim$A1A2.x == "AT" | outcome_MCP_QC_bim$A1A2.x == "CG"]
remove_rsids_MCP <- c(remove_rsids_MCP, outcome_MCP_QC_bim$id[outcome_MCP_QC_bim$A1A2.x != outcome_MCP_QC_bim$A1A2.y])

remove_rsids_CD40 <- exposure_CD40_QC_bim$id[exposure_CD40_QC_bim$A1A2.x == "AT" | exposure_CD40_QC_bim$A1A2.x == "CG"]
remove_rsids_CD40 <- c(remove_rsids_CD40, exposure_CD40_QC_bim$id[exposure_CD40_QC_bim$A1A2.x != exposure_CD40_QC_bim$A1A2.y])

## Flip alleles, beta, zscore and freq for misaligned SNPs
outcome_MDD_QC_harmonised <- outcome_MDD_QC[(outcome_MDD_QC$id %in% bim_file$id) & (!outcome_MDD_QC$id %in% remove_rsids_MDD),]
outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "beta"] <- outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "beta"] * -1
outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "zscore"] <- outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "zscore"] * -1
outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "freq"] <- 1 - outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "freq"]
## Duplicate A1 column so information is not overwritten
outcome_MDD_QC_harmonised$A1_spare <- outcome_MDD_QC_harmonised$A1
outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "A1"] <- outcome_MDD_QC_harmonised$A2[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD]
outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "A2"] <- outcome_MDD_QC_harmonised$A1_spare[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD]

outcome_MCP_QC_harmonised <- outcome_MCP_QC[(outcome_MCP_QC$id %in% bim_file$id) & (!outcome_MCP_QC$id %in% remove_rsids_MCP),]
outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "beta"] <- outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "beta"] * -1
outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "zscore"] <- outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "zscore"] * -1
outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "freq"] <- 1 - outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "freq"]
## Duplicate A1 column so information is not overwritten
outcome_MCP_QC_harmonised$A1_spare <- outcome_MCP_QC_harmonised$A1
outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "A1"] <- outcome_MCP_QC_harmonised$A2[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP]
outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "A2"] <- outcome_MCP_QC_harmonised$A1_spare[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP]

exposure_CD40_QC_harmonised <- exposure_sCD40_QC[(exposure_sCD40_QC$id %in% bim_file$id) & (!exposure_sCD40_QC$id %in% remove_rsids_CD40),]
exposure_CD40_QC_harmonised[exposure_CD40_QC_harmonised$id %in% flipped_rsids_CD40, "beta"] <- exposure_CD40_QC_harmonised[exposure_CD40_QC_harmonised$id %in% flipped_rsids_CD40, "beta"] * -1
exposure_CD40_QC_harmonised[exposure_CD40_QC_harmonised$id %in% flipped_rsids_CD40, "zscore"] <- exposure_CD40_QC_harmonised[exposure_CD40_QC_harmonised$id %in% flipped_rsids_CD40, "zscore"] * -1
exposure_CD40_QC_harmonised[exposure_CD40_QC_harmonised$id %in% flipped_rsids_CD40, "freq"] <- 1 - exposure_CD40_QC_harmonised[exposure_CD40_QC_harmonised$id %in% flipped_rsids_CD40, "freq"]
## Duplicate A1 column so information is not overwritten
exposure_CD40_QC_harmonised$A1_spare <- exposure_CD40_QC_harmonised$A1
exposure_CD40_QC_harmonised[exposure_CD40_QC_harmonised$id %in% flipped_rsids_CD40, "A1"] <- exposure_CD40_QC_harmonised$A2[exposure_CD40_QC_harmonised$id %in% flipped_rsids_CD40]
exposure_CD40_QC_harmonised[exposure_CD40_QC_harmonised$id %in% flipped_rsids_CD40, "A2"] <- exposure_CD40_QC_harmonised$A1_spare[exposure_CD40_QC_harmonised$id %in% flipped_rsids_CD40]

## Double check harmonisation has worked
bim_sumstats_check <- full_join(bim_file[,c("id", "pos", "A1", "A2")], (outcome_MDD_QC_harmonised %>% 
                                                                          select("id", "pos", "A1", "A2") %>%
                                                                          rename(A1_MDD = A1,
                                                                                 A2_MDD = A2)), by = c("id", "pos"))

bim_sumstats_check <- full_join(bim_sumstats_check, (outcome_MCP_QC_harmonised %>% 
                                                       select("id", "pos", "A1", "A2") %>%
                                                       rename(A1_MCP = A1,
                                                              A2_MCP = A2)), by = c("id", "pos"))


bim_sumstats_check <- full_join(bim_sumstats_check, (exposure_CD40_QC_harmonised %>% 
                                                       select("id", "pos", "A1", "A2") %>%
                                                       rename(A1_CD40 = A1,
                                                              A2_CD40 = A2)), by = c("id", "pos"))


table(bim_sumstats_check$A1 == bim_sumstats_check$A1_MDD)
table(bim_sumstats_check$A1 == bim_sumstats_check$A1_MCP)
table(bim_sumstats_check$A1 == bim_sumstats_check$A1_CD40)


# ## Load in BED reference file 
# x <- read_reference_bed("~/Desktop/PhD/resources/LD_reference/eur_CD40.bed", verbose = FALSE)
# 
# ## Harmonise MDD sumstats
# hsumstats_MDD <- harmonize_sumstats(outcome_MDD_QC, x, check_strand_flip = TRUE)
# temp <- left_join(hsumstats_MDD, outcome_MDD_QC, by = "id")
# flipped_rsids_MDD <- temp$id[temp$A1.x != temp$A1.y]
# 
# ## reverse beta and frequency for flipped alleles
# outcome_MDD_QC_harmonised <- outcome_MDD_QC[outcome_MDD_QC$id %in% hsumstats_MDD$id,]
# outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "beta"] <- outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "beta"] * -1
# outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "zscore"] <- outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "zscore"] * -1
# outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "freq"] <- 1 - outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "freq"]
# ## Duplicate A1 column so information is not overwritten
# outcome_MDD_QC_harmonised$A1_spare <- outcome_MDD_QC_harmonised$A1
# outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "A1"] <- outcome_MDD_QC_harmonised$A2[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD]
# outcome_MDD_QC_harmonised[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD, "A2"] <- outcome_MDD_QC_harmonised$A1_spare[outcome_MDD_QC_harmonised$id %in% flipped_rsids_MDD]
# 
# ## Harmonise MCP sumstats
# hsumstats_MCP<- harmonize_sumstats(outcome_MCP_QC, x, check_strand_flip = TRUE)
# temp <- inner_join(hsumstats_MCP, outcome_MCP_QC, by = "id")
# flipped_rsids_MCP <- temp$id[temp$A1.x != temp$A1.y]
# 
# ## reverse beta and frequency for flipped alleles
# outcome_MCP_QC_harmonised <- outcome_MCP_QC[outcome_MCP_QC$id %in% hsumstats_MCP$id,]
# outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "beta"] <- outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "beta"] * -1
# outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "zscore"] <- outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "zscore"] * -1
# outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "freq"] <- 1 - outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "freq"]
# ## Duplicate A1 column so information is not overwritten
# outcome_MCP_QC_harmonised$A1_spare <- outcome_MCP_QC_harmonised$A1
# outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "A1"] <- outcome_MCP_QC_harmonised$A2[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP]
# outcome_MCP_QC_harmonised[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP, "A2"] <- outcome_MCP_QC_harmonised$A1_spare[outcome_MCP_QC_harmonised$id %in% flipped_rsids_MCP]
# 
# ## Harmonise sCD40 sumstats
# hsumstats_sCD40 <- harmonize_sumstats(exposure_sCD40_QC, x, check_strand_flip = FALSE)
# temp <- left_join(hsumstats_sCD40, exposure_sCD40_QC, by = "id")
# flipped_rsids_sCD40 <- temp$id[temp$A1.x != temp$A1.y]
#  
# ## reverse beta and frequency for flipped alleles
# exposure_sCD40_QC_harmonised <- exposure_sCD40_QC[exposure_sCD40_QC$id %in% hsumstats_sCD40$id,]
# exposure_sCD40_QC_harmonised[exposure_sCD40_QC_harmonised$id %in% flipped_rsids_sCD40, "beta"] <- exposure_sCD40_QC_harmonised[exposure_sCD40_QC_harmonised$id %in% flipped_rsids_sCD40, "beta"] * -1
# exposure_sCD40_QC_harmonised[exposure_sCD40_QC_harmonised$id %in% flipped_rsids_sCD40, "zscore"] <- exposure_sCD40_QC_harmonised[exposure_sCD40_QC_harmonised$id %in% flipped_rsids_sCD40, "zscore"] * -1
# exposure_sCD40_QC_harmonised[exposure_sCD40_QC_harmonised$id %in% flipped_rsids_sCD40, "freq"] <- 1 - exposure_sCD40_QC_harmonised[exposure_sCD40_QC_harmonised$id %in% flipped_rsids_sCD40, "freq"]
# ## Duplicate A1 column so information is not overwritten
# exposure_sCD40_QC_harmonised$A1_spare <- exposure_sCD40_QC_harmonised$A1
# exposure_sCD40_QC_harmonised[exposure_sCD40_QC_harmonised$id %in% flipped_rsids_sCD40, "A1"] <- exposure_sCD40_QC_harmonised$A2[exposure_sCD40_QC$id %in% flipped_rsids_sCD40]
# exposure_sCD40_QC_harmonised[exposure_sCD40_QC_harmonised$id %in% flipped_rsids_sCD40, "A2"] <- exposure_sCD40_QC_harmonised$A1_spare[exposure_sCD40_QC$id %in% flipped_rsids_sCD40]

# ## Remove ambiguous SNPs
# ambiguous_SNPs <- exposure_sCD40_QC_harmonised$id[exposure_sCD40_QC_harmonised$A1 == "T" & exposure_sCD40_QC_harmonised$A2 == "A"] 
# ambiguous_SNPs <- c(ambiguous_SNPs, exposure_sCD40_QC_harmonised$id[exposure_sCD40_QC_harmonised$A1 == "A" & exposure_sCD40_QC_harmonised$A2 == "T"])
# ambiguous_SNPs <- c(ambiguous_SNPs, exposure_sCD40_QC_harmonised$id[exposure_sCD40_QC_harmonised$A1 == "C" & exposure_sCD40_QC_harmonised$A2 == "G"])
# ambiguous_SNPs <- c(ambiguous_SNPs, exposure_sCD40_QC_harmonised$id[exposure_sCD40_QC_harmonised$A1 == "G" & exposure_sCD40_QC_harmonised$A2 == "C"])
# 
# exposure_sCD40_QC_harmonised_no_ambiguous <- exposure_sCD40_QC_harmonised[(!exposure_sCD40_QC_harmonised$id %in% ambiguous_SNPs)]

## Save harmonised sumstats
write.table(outcome_MDD_QC_harmonised, '~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/MDD_harmonised.tsv', sep = '\t', row.names = F, quote = F)
write.table(outcome_MCP_QC_harmonised, '~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/MCP_harmonised.tsv', sep = '\t', row.names = F, quote = F)
write.table(exposure_CD40_QC_harmonised, '~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/sCD40_harmonised.tsv', sep = '\t', row.names = F, quote = F)
