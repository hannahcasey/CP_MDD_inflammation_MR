MDD <- read.table('~/Desktop/PhD/resources/summary_statistics/processed_data/MDD.txt', header = T, fill = T)
MDD_old <- read.table('~/Desktop/PhD/resources/summary_statistics/raw_data/MDD_old.txt', header = T, fill = T)

MCP <- read.table("/Users/hannahcasey/Desktop/PhD/resources/summary_statistics/raw_data/ChronicPain.stats", header = T)
CD40 <- read.table("~/Desktop/PhD/resources/summary_statistics/raw_data/olink_CVD/TNFRSF5.txt", header = TRUE)
  
## Get cis region of exposure protein
## Chromosome 20
bp_range = (46118316 - 500000):(46129863  + 500000)

## Get SNPs in cis region in exposure and outcome GWAS
# Get SNPs within 1Mb of independently associated SNPs (p > 5x10-8)
CD40_snps <- sCD40$hm_rsid[(sCD40$hm_chrom == 20 & sCD40$hm_pos %in% bp_range)]

## Get SNPs in in both exposure and outcome GWAS
MDD_snps <- CD40_snps[CD40_snps %in% MDD$SNP]
MCP_snps <- CD40_snps[CD40_snps %in% MCP$SNP]

## Combine snp lists
TNFRSF5_common_snps <- unique(c(MDD_snps, MCP_snps))

## write SNP list to calculate LD matrix
write.table(TNFRSF5_common_snps, "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/TNFRSF5_common_snps.txt", row.names = F, quote = F, col.names = F)
