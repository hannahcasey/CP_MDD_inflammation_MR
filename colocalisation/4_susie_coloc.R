library(coloc)
library(dplyr)
#install.packages("devtools")
library(devtools)
#install_github("jrs95/gassocplot")
library(gassocplot)
#install.packages("LDlinkR")
library(LDlinkR)
library(ieugwasr)
#install.packages("Rfast")
library(Rfast)
library(patchwork)
library(ggplot2)
#install.packages("susieR")
library(susieR)
set.seed(151097)
library(LDlinkR)
library(ggplotify)

## Read in harmonized summary statistics - limited to 500kb +/- CD40 transcipt region
MCP <- read.table('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/MCP_harmonised.tsv', header = T)
MDD <- read.table('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/MDD_harmonised.tsv', header = T)
CD40 <- read.table('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/sCD40_harmonised.tsv', header = T)


## Load LD matrix and SNP list
LD <- read.table('/Users/hannahcasey/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/LD_matrix.ld')
LD_snps <- read.table('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/LD_matrix.snplist')

## Find duplicates in SNP list
LD_no_dup <- LD[-which(duplicated(LD_snps$V1)), -which(duplicated(LD_snps$V1))]
LD_snps_no_dup <- LD_snps$V1[-which(duplicated(LD_snps$V1))]

colnames(LD_no_dup) <- LD_snps_no_dup
rownames(LD_no_dup) <- LD_snps_no_dup

## Remove rows with all NaN values from LD matrix
LD_no_dup <- LD_no_dup[-which(is.nan(as.matrix(LD_no_dup$rs143028411))), -which(is.nan(as.matrix(LD_no_dup$rs143028411)))]

# ## Remove SNPs with no LD info from GWAS sumstats
MCP <- MCP[MCP$id %in% LD_snps_no_dup,]
MDD <- MDD[MDD$id %in% LD_snps_no_dup,]
CD40 <- CD40[CD40$id %in% LD_snps_no_dup,]

## Get LD matrices for each GWAS
LD_MCP <- LD_no_dup[MCP$id[MCP$id %in% colnames(LD_no_dup)], MCP$id[MCP$id %in% colnames(LD_no_dup)]]
LD_MDD <- LD_no_dup[MDD$id[MDD$id %in% colnames(LD_no_dup)], MDD$id[MDD$id %in% colnames(LD_no_dup)]]
LD_CD40 <- LD_no_dup[CD40$id[CD40$id %in% colnames(LD_no_dup)], CD40$id[CD40$id %in% colnames(LD_no_dup)]]

## Plot CD40 region 
CD40_small <- CD40[CD40$id %in% colnames(LD_CD40),]
CD40_small_plot <-  CD40_small %>% dplyr::select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)
assoc_plot(CD40_small_plot, corr = LD_CD40)

MDD_small <- MDD[MDD$id %in% colnames(LD_MDD),]
assoc_plot((MDD_small %>% select(id, chr, pos, zscore) %>% rename(marker = id, z = zscore)), corr = LD_MDD)

MCP_small <- MCP[MCP$id %in% colnames(LD_MCP),]
assoc_plot((MCP_small %>% select(id, chr, pos, zscore) %>% rename(marker = id, z = zscore)), corr = LD_MCP)

## get MAF
CD40_small$MAF <- as.numeric(CD40_small$freq)
#CD40_small$MAF <- ifelse(CD40_small$freq > 0.5, 1-CD40_small$freq, CD40_small$freq)

MDD_small$MAF <- as.numeric(MDD_small$freq)
#MDD_small$MAF <- ifelse(MDD_small$freq > 0.5, 1-MDD_small$freq, MDD_small$freq)

MCP_small$MAF <- as.numeric(MCP_small$freq)
#MCP_small$MAF <- ifelse(MCP_small$freq > 0.5, 1-MCP_small$freq, MCP_small$freq)

## make list
exposure_list <- list("pvalues" = CD40_small$pvalue, "beta" = CD40_small$beta, "varbeta" = CD40_small$varbeta, "N" = CD40_small$n, "type" = "quant", "MAF" = CD40_small$MAF, "LD" = as.matrix(LD_CD40), "snp" = CD40_small$id, "position" = CD40_small$pos)
outcome_MDD_list <- list("pvalues" = MDD_small$pvalue, "beta" = MDD_small$beta, "varbeta" = MDD_small$varbeta, "N" = MDD_small$n, "type" = "cc", "s" = mean(MDD_small$s), "MAF" = MDD_small$MAF, "LD" = as.matrix(LD_MDD), "snp" = MDD_small$id, "position" = MDD_small$pos, "sdY" = 1)
outcome_MCP_list <- list("pvalues" = MCP_small$pvalue, "beta" = MCP_small$beta, "varbeta" = MCP_small$varbeta, "N" = 440780, "type" = "quant", "MAF" = MCP_small$MAF, "LD" = as.matrix(LD_MCP), "snp" = MCP_small$id, "position" = MCP_small$pos)

## Check allele mismatch
condz_CD40 = kriging_rss(z = CD40_small$zscore, R = as.matrix(LD_CD40), n=mean(CD40_small$n))
condz_CD40$plot

## Remove SNPs with large difference between observed zscores and expected zscores (indication of LD panel mismatch)
remove_index <- which(abs(condz_CD40$conditional_dist$z_std_diff) > 2.5)
CD40_small_alligned <- CD40_small[-remove_index,]
LD_CD40_alligned <- LD_CD40[-remove_index, -remove_index]

## check again
condz_CD40_alligned = kriging_rss(z = CD40_small_alligned$zscore, R = as.matrix(LD_CD40_alligned), n=mean(CD40_small_alligned$n))
condz_CD40_alligned$plot

condz_MDD = kriging_rss(z = MDD_small$zscore, R = as.matrix(LD_MDD), n=mean(MDD_small$n))
condz_MDD$plot

condz_MCP = kriging_rss(z = MCP_small$zscore, R = as.matrix(LD_MCP), n=440780)
condz_MCP$plot

## Remake list
exposure_alligned_list <- list("pvalues" = CD40_small_alligned$pvalue, "beta" = CD40_small_alligned$beta, "varbeta" = CD40_small_alligned$varbeta, "N" = CD40_small_alligned$n, "type" = "quant", "MAF" = CD40_small_alligned$MAF, "LD" = as.matrix(LD_CD40_alligned), "snp" = CD40_small_alligned$id, "position" = CD40_small_alligned$pos)

## Check alignment - should have positive skew
check_alignment(exposure_alligned_list)
check_alignment(outcome_MDD_list)
check_alignment(outcome_MCP_list)

## standard coloc (single causal variant)
my.res <- coloc.abf(dataset1=exposure_alligned_list, dataset2=outcome_MDD_list)
sensitivity(my.res,"H4 > 0.9")

my.res <- coloc.abf(dataset1=exposure_alligned_list, dataset2=outcome_MCP_list)
sensitivity(my.res,"H4 > 0.9")

## finemap
finemap.signals(exposure_alligned_list,method="cond")
finemap.signals(outcome_MDD_list,method="cond")
finemap.signals(outcome_MCP_list,method="cond")

S_outcome_MDD_list <- runsusie(outcome_MDD_list, n = 3583863, max_causal = 100, max_iter = 1000, r2.prune = 0.2, coverage=0.01)
summary(S_outcome_MDD_list)

S_outcome_MCP_list <- runsusie(outcome_MCP_list, n = 440780, max_causal = 100, max_iter = 1000, coverage=0.001, r2.prune = 0.2)
summary(S_outcome_MCP_list)

S_exposure_list <- runsusie(exposure_alligned_list, n = 21758,  max_causal = 100, max_iter = 10000, r2.prune = 0.2, coverage=0.95, refine = TRUE)
summary(S_exposure_list)

## Plot each causal SNP
for (snp in (as.character(sapply(S_outcome_MCP_list[["sets"]]$cs,names)))){
  assoc_plot((MCP_small %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_MCP, top.marker=snp)
}

for (snp in (as.character(sapply(S_outcome_MDD_list[["sets"]]$cs,names)))){
  assoc_plot((MDD_small %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_MDD, top.marker=snp)
}

for (snp in (as.character(sapply(S_exposure_list[["sets"]]$cs,names)))){
   assoc_plot((CD40_small %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_CD40, top.marker="rs190760567")
}


## Save susie output
save(S_outcome_MDD_list, file = "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/S_outcome_MDD_list.Rdata")
save(S_outcome_MCP_list, file = "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/S_outcome_MCP_list.Rdata")
save(S_exposure_list,file = "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/S_exposure_list.Rdata")

# load("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/S_outcome_MDD_list.Rdata")
# load("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/S_outcome_MCP_list.Rdata")
# load("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/S_exposure_list.Rdata")


temp = coloc.susie(S_outcome_MCP_list,S_exposure_list)

## Colocalise CD40 and MCP
res_susie_CD40_MCP <- if(requireNamespace("susieR",quietly=TRUE)) {
  susie_MCP_res=coloc.susie(S_outcome_MCP_list,S_exposure_list)
  print(susie_MCP_res$summary)
}

for (i in 1:nrow(susie_MCP_res$summary)){
  
  if(requireNamespace("susieR",quietly=TRUE)) {
    sensitivity(susie_MCP_res,"H4 > 0.8",row=i,dataset1=outcome_MCP_list,dataset2=exposure_alligned_list)
  }
}


## Colocalise CD40 and MDD
res_susie_CD40_MDD <- if(requireNamespace("susieR",quietly=TRUE)) {
  susie_MDD_res=coloc.susie(S_outcome_MDD_list,S_exposure_list)
  print(susie_MDD_res$summary)
}

## Sensitivity analysis
for (i in 1:nrow(susie_MDD_res$summary)){
  
  if(requireNamespace("susieR",quietly=TRUE)) {
    sensitivity(susie_MDD_res,"H4 > 0.8",row=i,dataset1=outcome_MDD_list,dataset2=exposure_list)
  }
}

## Plot colocalised SNPs
for(i in 1:nrow(susie_MDD_res$summary)){
  
  snp1 <- susie_MDD_res$summary$hit1[i]
  snp2 <- susie_MDD_res$summary$hit2[i]
  
  p1 <- as.ggplot(assoc_plot((MDD_small %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_MDD, top.marker=snp1))
  p2 <- as.ggplot(assoc_plot((CD40_small %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_CD40, top.marker=snp2))
  
  assign(paste0("CD40_MDD_", i), p1/p2)
}

for(i in 1:nrow(susie_MCP_res$summary)){
  
  snp1 <- susie_MCP_res$summary$hit1[i]
  snp2 <- susie_MCP_res$summary$hit2[i]
  
  p1 <- as.ggplot(assoc_plot((MCP_small %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_MCP, top.marker=snp1))
  p2 <- as.ggplot(assoc_plot((CD40_small %>% select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_CD40, top.marker=snp2))
  
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


## plot LD matrix between causal SNPs
corrplot::corrplot(as.matrix(LD_CD40[unique(susie_MCP_res$summary$hit2), unique(susie_MCP_res$summary$hit2)]))

LD_CD40[unique(susie_MCP_res$summary$hit1), unique(susie_MCP_res$summary$hit1)]

