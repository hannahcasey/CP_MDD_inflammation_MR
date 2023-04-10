MDD <- read.table('~/Desktop/PhD/resources/summary_statistics/processed_data/MDD.txt', header = T, fill = T)
MCP <- read.table("/Users/hannahcasey/Desktop/PhD/resources/summary_statistics/raw_data/ChronicPain.stats", header = T)
CD40 <- read.table("/Users/hannahcasey/Desktop/PhD/resources/summary_statistics/raw_data/sCD40.txt", header = TRUE)


## limit FAM file to those of GBR ancestry
EUR <- read.csv("~/Desktop/PhD/resources/LD_reference/EUR_1000g_IDs.csv", header = F)
FAM_1000g <- read.table("~/Desktop/PhD/resources/LD_reference/all_phase3.fam", header = FALSE)
FAM_1000g_EUR <- FAM_1000g[FAM_1000g$V2 %in% EUR$V2,]

write.table(FAM_1000g, "~/Desktop/PhD/resources/LD_reference/all_phase3_original.fam", quote = FALSE, row.names = FALSE)
write.table(FAM_1000g_EUR, "~/Desktop/PhD/resources/LD_reference/all_phase3_eur.fam", quote = FALSE, row.names = FALSE)

bp_range_hg38 = (46118316 - 500000):(46129863  + 500000)
bp_range_hg19 = ( 44746910 - 500000):(44758502  + 500000)

## Get SNPs in cis region in exposure and outcome GWAS
CD40_snps <- CD40$hm_rsid[(CD40$hm_chrom == 20 & CD40$hm_pos %in% bp_range_hg38)]
MDD_snps <- MDD$ID[(MDD$CHROM == 20 & MDD$POS %in% bp_range_hg19)]
MCP_snps <- MCP$SNP[(MCP$CHR == 20 & MCP$BP %in% bp_range_hg19)]

## Combine snp lists
CD40_common_snps <- unique(c(CD40_snps, MDD_snps, MCP_snps))

## write SNP list to calculate LD matrix
write.table(CD40_common_snps, "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/CD40_common_snps.txt", row.names = F, quote = F, col.names = F)
