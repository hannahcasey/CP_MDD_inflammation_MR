library(coloc)
library(dplyr)
#install.packages("devtools")
library(devtools)
#install_github("jrs95/gassocplot")
library(gassocplot)
library(patchwork)
library(ggplot2)
library(susieR)
set.seed(151097)
library(ggplotify)

## Read in harmonized summary statistics - limited to 500kb +/- CD40 transcript region
MCP <- read.table('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/MCP_CD40_harmonised.tsv', header = T)
MDD <- read.table('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/MDD_CD40_harmonised.tsv', header = T)
CD40 <- read.table('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/sCD40_CD40_harmonised.tsv', header = T)

## Replace smallest p-values with smallest number possible in R
CD40$pvalue[CD40$pvalue == 0] <- 2.225074e-308

## Load LD matrix and SNP list
LD <- read.table('/Users/hannahcasey/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/UKB_CD40_LD.ld')
LD_snps <- read.table('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/UKB_CD40_LD.snplist')

## Find duplicates in SNP list
LD_no_dup <- LD[-which(duplicated(LD_snps$V1)), -which(duplicated(LD_snps$V1))]
LD_snps_no_dup <- LD_snps$V1[-which(duplicated(LD_snps$V1))]

colnames(LD_no_dup) <- LD_snps_no_dup
rownames(LD_no_dup) <- LD_snps_no_dup

## Remove SNPs with missing NAs from LD
LD_no_dup_complete <- LD_no_dup[LD_snps_no_dup[complete.cases(LD_no_dup)], LD_snps_no_dup[complete.cases(LD_no_dup)]]

## Remove SNPs with no LD info from GWAS sumstats
MCP <- MCP[MCP$id %in% names(LD_no_dup_complete),]
MDD <- MDD[MDD$id %in% names(LD_no_dup_complete),]
CD40 <- CD40[CD40$id %in% names(LD_no_dup_complete),]

## Only keep SNPs common to both exposure (sCD40) and outcomes (MCP and MDD)
MCP <- MCP[MCP$id %in% CD40$id, ]
MDD <- MDD[MDD$id %in% CD40$id, ]
CD40_MCP <- CD40[CD40$id %in% MCP$id, ]
CD40_MDD <- CD40[CD40$id %in% MDD$id, ]

## get MAF
MDD$MAF <- as.numeric(MDD$freq)
MDD$MAF[MDD$MAF > 0.5] <- 1 - MDD$MAF[MDD$MAF > 0.5]

MCP$MAF <- as.numeric(MCP$freq)
MCP$MAF[MCP$MAF > 0.5] <- 1 - MCP$MAF[MCP$MAF > 0.5]

CD40_MDD$MAF <- as.numeric(CD40_MDD$freq)
CD40_MDD$MAF[CD40_MDD$MAF > 0.5] <- 1 - CD40_MDD$MAF[CD40_MDD$MAF > 0.5]

CD40_MCP$MAF <- as.numeric(CD40_MCP$freq)
CD40_MCP$MAF[CD40_MCP$MAF > 0.5] <- 1 - CD40_MCP$MAF[CD40_MCP$MAF > 0.5]

## Get LD matrices for each sumstat
LD_MDD <- LD_no_dup_complete[MDD$id,MDD$id]
LD_MCP <- LD_no_dup_complete[MCP$id,MCP$id]
LD_CD40_MDD <- LD_no_dup_complete[CD40_MDD$id,CD40_MDD$id]
LD_CD40_MCP <- LD_no_dup_complete[CD40_MCP$id,CD40_MCP$id]

LD_CD40 <- LD_no_dup_complete[CD40$id,CD40$id]

## Reorder sumstats and LD by position
MDD <- MDD[order(MDD$pos),]
MCP <- MCP[order(MCP$pos),]
CD40_MDD <- CD40_MDD[order(CD40_MDD$pos),]
CD40_MCP <- CD40_MCP[order(CD40_MCP$pos),]

CD40 <- CD40[order(CD40$pos),]

LD_MDD <- LD_MDD[MDD$id, MDD$id]
LD_MCP <- LD_MCP[MCP$id, MCP$id]
LD_CD40_MDD <- LD_CD40_MDD[CD40_MDD$id, CD40_MDD$id]
LD_CD40_MCP <- LD_CD40_MCP[CD40_MCP$id, CD40_MCP$id]

LD_CD40 <- LD_CD40[CD40$id, CD40$id]

## Check for allele flip issues ----
RSS_MDD <- kriging_rss(MDD$zscore, LD_MDD, n = median(MDD$n))
RSS_MDD$plot

RSS_MCP <- kriging_rss(MCP$zscore, LD_MCP, n = 440780)
RSS_MCP$plot

RSS_CD40_MDD <- kriging_rss(CD40_MDD$zscore, LD_CD40_MDD, n = median(CD40_MDD$n))
RSS_CD40_MDD$plot

RSS_CD40_MCP <- kriging_rss(CD40_MCP$zscore, LD_CD40_MCP, n = median(CD40_MCP$n))
RSS_CD40_MCP$plot

# flipped_SNP_index <- which(with(RSS_CD40$conditional_dist, logLR>2 & abs(z)>2))
# CD40 <- CD40[-flipped_SNP_index,]
# LD_CD40 <- LD_CD40[-flipped_SNP_index, -flipped_SNP_index]
# RSS_CD40 <- kriging_rss(CD40$zscore, LD_CD40, n = median(CD40$n))
# flipped_SNP_index  <- which(with(RSS_CD40$conditional_dist, logLR>2 & abs(z)>2))

## Make list to run coloc
MCP_list <- list("pvalues" = MCP$pvalue, "beta" = MCP$beta, "varbeta" = MCP$varbeta, "N" = 440780, "type" = "quant", "MAF" = MCP$MAF, "LD" = as.matrix(LD_MCP), "snp" = MCP$id, "position" = MCP$pos)
MDD_list <- list("pvalues" = MDD$pvalue, "beta" = MDD$beta, "varbeta" = MDD$varbeta, "N" = MDD$n, "type" = "cc", "s" = median(MDD$s), "MAF" = MDD$MAF, "snp" = MDD$id, "position" = MDD$pos, "LD" = as.matrix(LD_MDD))
CD40_MDD_list <- list("pvalues" = CD40_MDD$pvalue, "beta" = CD40_MDD$beta, "varbeta" = CD40_MDD$varbeta, "N" = CD40_MDD$n, "type" = "quant", "MAF" = CD40_MDD$MAF, "snp" = CD40_MDD$id, "position" = CD40_MDD$pos, "LD" = as.matrix(LD_CD40_MDD))
CD40_MCP_list <- list("pvalues" = CD40_MCP$pvalue, "beta" = CD40_MCP$beta,  "varbeta" = CD40_MCP$varbeta, "N" = CD40_MCP$n, "type" = "quant", "MAF" = CD40_MCP$MAF, "snp" = CD40_MCP$id, "position" = CD40_MCP$pos, "LD" = as.matrix(LD_CD40_MCP))

# MDD_list <- list("pvalues" = MDD$pvalue, "N" = MDD$n, "type" = "cc", "s" = median(MDD$s), "MAF" = MDD$MAF, "snp" = MDD$id, "position" = MDD$pos, "LD" = as.matrix(LD_MDD))
# CD40_MDD_list <- list("pvalues" = CD40_MDD$pvalue, "N" = CD40_MDD$n, "type" = "quant", "MAF" = CD40_MDD$MAF, "snp" = CD40_MDD$id, "position" = CD40_MDD$pos, "LD" = as.matrix(LD_CD40_MDD))
# 
# MDD_list <- list("pvalues" = MDD$pvalue, "beta" = MDD$beta, "varbeta" = MDD$varbeta, "type" = "cc", "LD" = as.matrix(LD_MDD), "snp" = MDD$id, "position" = MDD$pos)
# CD40_MDD_list <- list("pvalues" = CD40_MDD$pvalue, "beta" = CD40_MDD$beta, "varbeta" = CD40_MDD$varbeta, "sdY" = 1, "N" = CD40_MDD$n, "type" = "quant", "snp" = CD40_MDD$id, "position" = CD40_MDD$pos, "LD" = as.matrix(LD_CD40_MDD))

## Run standard coloc and sensitivity analysis
CD40_MDD_res <- coloc.abf(dataset1=CD40_MDD_list, dataset2=MDD_list)
## Change PP.H4 calculated as 0 to lowest number possible, otherwise no color is assigned for plotting
CD40_MDD_res$results$SNP.PP.H4[CD40_MDD_res$results$SNP.PP.H4 == 0] <- 2.225074e-308
coloc::sensitivity(CD40_MDD_res,"H4 > 0.9")

CD40_MCP_res <- coloc.abf(dataset1=CD40_MCP_list, dataset2=MCP_list)
## Change PP.H4 calculated as 0 to lowest number possible, otherwise no color is assigned for plotting
CD40_MCP_res$results$SNP.PP.H4[CD40_MCP_res$results$SNP.PP.H4 == 0] <- 2.225074e-308
coloc::sensitivity(CD40_MCP_res,"H4 > 0.9")

# ## Plot coloc using with LD
top_snp <- CD40_MDD_res$results$snp[which(CD40_MDD_res$results$SNP.PP.H4 == max(CD40_MDD_res$results$SNP.PP.H4))]
p1 <- as.ggplot(assoc_plot((MDD %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_MDD, top.marker= top_snp))
p2 <- as.ggplot(assoc_plot((CD40_MDD %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_CD40_MDD, top.marker=top_snp))
p1/p2

top_snp <- CD40_MCP_res$results$snp[which(CD40_MCP_res$results$SNP.PP.H4 == max(CD40_MCP_res$results$SNP.PP.H4))]
p1 <- as.ggplot(assoc_plot((MCP %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_MCP, top.marker= top_snp))
p2 <- as.ggplot(assoc_plot((CD40_MCP %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_CD40_MCP, top.marker=top_snp))
p1/p2

## Check alignment - should have positive skew (does not work with varbeta removed)
check_alignment(CD40_MCP_list)
check_alignment(CD40_MDD_list)
check_alignment(MDD_list)
check_alignment(MCP_list)

## Carry out susie finemapping
MDD_susie <- runsusie(MDD_list, n = 3583863, max_causal = 100, max_iter = 1000, coverage=0.95)
summary(MDD_susie)

MCP_susie <- runsusie(MCP_list, n = 440780, max_causal = 100, max_iter = 1000, coverage=0.1)
summary(MCP_susie)

CD40_MDD_susie <- runsusie(CD40_MDD_list, n = 21758,  max_causal = 100, max_iter = 1000, coverage=0.95)
summary(CD40_MDD_susie)

CD40_MCP_susie <- runsusie(CD40_MCP_list, n = 21758,  max_causal = 100, max_iter = 1000, coverage=0.95)
summary(CD40_MCP_susie)

## Colocalise CD40 and MCP
res_susie_CD40_MCP <- if(requireNamespace("susieR",quietly=TRUE)) {
  susie_MCP_res=coloc.susie(CD40_MCP_susie, MCP_susie)
  print(susie_MCP_res$summary)
}

for (i in 1:nrow(susie_MCP_res$summary)){
  
  susie_MCP_res$results[[paste0("SNP.PP.H4.row", i)]][susie_MCP_res$results[[paste0("SNP.PP.H4.row", i)]] == 0] <- 2.225074e-308

  if(requireNamespace("susieR",quietly=TRUE)) {
    sensitivity(susie_MCP_res,"H4 > 0.9",row=i,dataset1=MCP_list,dataset2=CD40_MCP_list)
  }
}

## Colocalise CD40 and MDD
res_susie_CD40_MDD <- if(requireNamespace("susieR",quietly=TRUE)) {
  susie_MDD_res=coloc.susie(MDD_susie, CD40_MDD_susie)
  print(susie_MDD_res$summary)
}

## Sensitivity analysis
for (i in 1:nrow(susie_MDD_res$summary)){
  
  susie_MDD_res$results[[paste0("SNP.PP.H4.row", i)]][susie_MDD_res$results[[paste0("SNP.PP.H4.row", i)]] == 0] <- 2.225074e-308
  
  if(requireNamespace("susieR",quietly=TRUE)) {
    sensitivity(susie_MDD_res,"H4 > 0.9",row=i,dataset1=MDD_list,dataset2=CD40_MDD_list)
  }
}

## Plot colocalised SNPs
for(i in 1:nrow(susie_MDD_res$summary)){

  snp1 <- susie_MDD_res$summary$hit1[i]
  snp2 <- susie_MDD_res$summary$hit2[i]

  p1 <- as.ggplot(assoc_plot((MDD %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_CD40_MDD, top.marker=snp1))
  p2 <- as.ggplot(assoc_plot((CD40_MDD %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_CD40_MDD, top.marker=snp2))

  assign(paste0("CD40_MDD_", i), p1/p2)
}

for(i in 1:nrow(susie_MCP_res$summary)){

  snp1 <- susie_MCP_res$summary$hit1[i]
  snp2 <- susie_MCP_res$summary$hit2[i]

  p1 <- as.ggplot(assoc_plot((MCP %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_MCP, top.marker=snp2))
  p2 <- as.ggplot(assoc_plot((CD40_MCP %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_CD40_MCP, top.marker=snp1))

  assign(paste0("CD40_MCP_", i), p1/p2)
}

## Save results
save(susie_MDD_res, file = "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/colocalization/susie_MDD_res.Rdata")
write.csv(susie_MDD_res$summary, file = "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/colocalization/susie_MDD_res.csv")

save(susie_MCP_res, file = "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/colocalization/susie_MCP_res.Rdata")
write.csv(susie_MCP_res$summary, file = "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/colocalization/susie_MCP_res.csv")

plots <- ls(patter = "CD40_MDD_|CD40_MCP_")

for (i in plots){
  png(paste0("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/colocalization/locuszoom/",i, ".jpg"), height = 1000)
  print(get(i))
  dev.off()
}
