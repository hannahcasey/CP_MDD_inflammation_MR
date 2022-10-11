library(coloc)
library(dplyr)
library(curl)
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
#install.packages("ggplotify")
library(ggplotify)
set.seed(151097)

## Read in harmonised summary statistics - limited to 500kb +/- CD40 transcipt region
MCP <- read.table('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/MCP_harmonised.tsv', header = T)
MDD <- read.table('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/MDD_harmonised.tsv', header = T)
CD40 <- read.table('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/CD40_harmonised.tsv', header = T)

## Load LD matric and SNP list
LD <- read.table('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/LD_matrix.ld')
LD_snps <- read.table('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/LD_matrix.snplist')

## find duplicates in SNP list
LD_no_dup <- LD[-which(duplicated(LD_snps$V1)), -which(duplicated(LD_snps$V1))]
LD_snps_no_dup <- LD_snps$V1[-which(duplicated(LD_snps$V1))]

colnames(LD_no_dup) <- LD_snps_no_dup
rownames(LD_no_dup) <- LD_snps_no_dup

LD_MCP <- LD_no_dup[colnames(LD_no_dup) %in% MCP$id,]
LD_MCP <- LD_MCP[,colnames(LD_no_dup) %in% MCP$id]

LD_MDD <- LD_no_dup[colnames(LD_no_dup) %in% MDD$id,]
LD_MDD <- LD_MDD[,colnames(LD_no_dup) %in% MDD$id]

LD_CD40 <- LD_no_dup[colnames(LD_no_dup) %in% CD40$id,]
LD_CD40 <- LD_CD40[,colnames(LD_no_dup) %in% CD40$id]

## Plot CD40 region 
CD40_small <- CD40[CD40$id %in% colnames(LD_CD40),]
assoc_plot((CD40_small %>% dplyr::select(id, chr, pos, zscore) %>% dplyr::rename(marker = id, z = zscore)), corr = LD_CD40)

MDD_small <- MDD[MDD$id %in% colnames(LD_MDD),]
assoc_plot((MDD_small %>% select(id, chr, pos, zscore) %>% rename(marker = id, z = zscore)), corr = LD_MDD)

MCP_no_dup[MCP_ %in% colnames(LD_MCP),]

MCP_small <- MCP[MCP$id %in% colnames(LD_MCP),]
assoc_plot((MCP_small %>% select(id, chr, pos, zscore) %>% rename(marker = id, z = zscore)), corr = LD_MCP)

## Evident from locus that mutlliple independent signals should be considered

## get MAF
CD40_small$freq <- as.numeric(CD40_small$freq)
CD40_small$MAF <- ifelse(CD40_small$freq > 0.5, 1-CD40_small$freq, CD40_small$freq)

MDD_small$freq <- as.numeric(MDD_small$freq)
MDD_small$MAF <- ifelse(MDD_small$freq > 0.5, 1-MDD_small$freq, MDD_small$freq)

MCP_small$freq <- as.numeric(MCP_small$freq)
MCP_small$MAF <- ifelse(MCP_small$freq > 0.5, 1-MCP_small$freq, MCP_small$freq)

## remove rows with NA
CD40_small <- na.omit(CD40_small)
MDD_small <- na.omit(MDD_small)
MCP_small <- na.omit(MCP_small)

## make list
exposure_list <- list("pvalues" = CD40_small$pvalue, "beta" = CD40_small$beta, "varbeta" = CD40_small$varbeta, "N" = 21758, "type" = "quant", "MAF" = CD40_small$MAF, "LD" = as.matrix(LD_CD40), "snp" = CD40_small$id, "position" = CD40_small$pos, "sdY" = 1)
outcome_MDD_list <- list("pvalues" = MDD_small$pvalue, "beta" = MDD_small$beta, "varbeta" = MDD_small$varbeta, "N" = 3583863, "type" = "quant", "MAF" = MDD_small$MAF, "LD" = as.matrix(LD_MDD), "snp" = MDD_small$id, "position" = MDD_small$pos, "sdY" = 1)
outcome_MCP_list <- list("pvalues" = MCP_small$pvalue, "beta" = MCP_small$beta, "varbeta" = MCP_small$varbeta, "N" = 440780, "type" = "quant", "MAF" = MCP_small$MAF, "LD" = as.matrix(LD_MCP), "snp" = MCP_small$id, "position" = MCP_small$pos, "sdY" = 1)

str(exposure_list)
str(outcome_MDD_list)
str(outcome_MCP_list)

check_dataset(exposure_list)
check_dataset(outcome_MDD_list)
check_dataset(outcome_MCP_list)

plot_dataset(exposure_list)
plot_dataset(outcome_MDD_list)
plot_dataset(outcome_MCP_list)


## standard coloc (single causal variant)
my.res <- coloc.abf(dataset1=exposure_list, dataset2=outcome_MDD_list)
sensitivity(my.res,"H4 > 0.9")

my.res <- coloc.abf(dataset1=exposure_list, dataset2=outcome_MCP_list)
sensitivity(my.res,"H4 > 0.9")


## finemap
finemap.signals(exposure_list,method="cond")
finemap.signals(outcome_MDD_list,method="cond")
finemap.signals(outcome_MCP_list,method="cond")

check_dataset(exposure_list,req="LD")
check_dataset(outcome_MDD_list,req="LD")
check_dataset(outcome_MCP_list,req="LD")

S_outcome_MDD_list <- runsusie(outcome_MDD_list, n = 3583863, max_causal = 10, max_iter = 1000)
summary(S_outcome_MDD_list)

S_outcome_MCP_list <- runsusie(outcome_MCP_list, n = 440780, max_causal = 10, max_iter = 1000)
summary(S_outcome_MCP_list)

S_exposure_list <- runsusie(exposure_list, n = 21758, max_causal = 10, max_iter = 1000)
summary(S_exposure_list)

 
## Save susie output
save(S_outcome_MDD_list, file = "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/S_outcome_MDD_list.Rdata")
save(S_outcome_MCP_list, file = "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/S_outcome_MCP_list.Rdata")
save(S_exposure_list,file = "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/S_exposure_list.Rdata")

# load("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/S_outcome_MDD_list.Rdata")
# load("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/S_outcome_MCP_list.Rdata")
# load("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/S_exposure_list.Rdata")

## Colocalise CD40 and MCP

res_susie_CD40_MCP <- if(requireNamespace("susieR",quietly=TRUE)) {
  susie_MCP_res=coloc.susie(S_outcome_MCP_list,S_exposure_list)
  print(susie_MCP_res$summary)
}

for (i in 1:nrow(susie_MCP_res$summary)){
  
  if(requireNamespace("susieR",quietly=TRUE)) {
    sensitivity(susie_MCP_res,"H4 > 0.9",row=i,dataset1=outcome_MCP_list,dataset2=exposure_list)
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
    sensitivity(susie_MDD_res,"H4 > 0.9",row=i,dataset1=outcome_MDD_list,dataset2=exposure_list)
  }
}

## Plot colocalisted SNPs
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

