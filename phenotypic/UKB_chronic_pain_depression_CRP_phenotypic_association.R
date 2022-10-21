## Load Packages
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
library(foreign)
#install.packages("ggplot2")
library(ggplot2)
library(dplyr)
#library(mosaic)
#install.packages("hrbrthemes")
library(hrbrthemes)
#library(jtools)
#install.packages("rstatix")
library(rstatix)
#install.packages("emmeans")
library(emmeans)
#install.packages("RNOmni")
library(RNOmni)
#install.packages("effects")
#install.packages("hrbrthemes")
library(hrbrthemes)
library(car) ## levenes test
library(tidyverse)

## Load preprocessed data from UKB_CP_MDD_inflam_prep.R
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
## Load in NMR metabolomics data
UKB_CP_MDD_inflam_covariates <- read.csv("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/resources/UKB_CP_MDD_inflam_covariates.csv", header = TRUE)

## Plot CRP in chronic pain groups/status and MDD status
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create temp df with no NAs in chronic pain status and CRP columns
UKB_CP_MDD_inflam_covariates_temp <- UKB_CP_MDD_inflam_covariates[!is.na(UKB_CP_MDD_inflam_covariates$chronic_pain_status) & !is.na(UKB_CP_MDD_inflam_covariates$CRP_log), ]

jpeg("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/phenotypic/plots/CRP_CP.jpeg")

## Plot CRP levels in chronic pain status
ggplot(UKB_CP_MDD_inflam_covariates_temp, aes(x=as.factor(chronic_pain_status), y= CRP_log, group = chronic_pain_status)) + 
  geom_boxplot(fill = c("#ffffff", "#2986cc"),alpha = 0.5) +
  theme_classic(base_size = 20) +
  xlab("Chronic Pain Group") + ylab("Serum CRP (log, mg/l)") +
  scale_x_discrete(labels=c("No Chronic Pain","Chronic Pain")) +
  ggtitle("")

dev.off()

## Create temp df with no NAs in chronic pain group info and CRP columns
UKB_CP_MDD_inflam_covariates_temp <- UKB_CP_MDD_inflam_covariates[!is.na(UKB_CP_MDD_inflam_covariates$chronic_pain_group) & !is.na(UKB_CP_MDD_inflam_covariates$CRP_log), ]

## Sort chronic pain groups prior to plotting:
UKB_CP_MDD_inflam_covariates_temp$chronic_pain_group <- factor(UKB_CP_MDD_inflam_covariates_temp$chronic_pain_group  , levels=c("No Sites" , "1-2 Sites", "3-4 Sites", "5-7 Sites", "Widespread"))

jpeg("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/phenotypic/plots/CRP_CP_sites.jpeg")

## Plot CRP levels in chronic pain site groups
ggplot(UKB_CP_MDD_inflam_covariates_temp, aes(x=chronic_pain_group, y= CRP_log, group = chronic_pain_group, fill = chronic_pain_group)) + 
  geom_boxplot(fill = c("#ffffff", "#B6DEFF", "#53B0FC", "#008CFE", "#006DC5"),alpha = 0.5) +
  theme_classic(base_size = 20) +
  labs(title="",
       x ="Chronic Pain Group", y = "Serum CRP (log, mg/l)")
dev.off()


## Plot CRP in recurrent MDD status groups
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create temp df with no NAs in chronic pain status and CRP columns
UKB_CP_MDD_inflam_covariates_temp <- UKB_CP_MDD_inflam_covariates[!is.na(UKB_CP_MDD_inflam_covariates$recurrent_depression) & !is.na(UKB_CP_MDD_inflam_covariates$CRP_log), ]

jpeg("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/phenotypic/plots/CRP_MDD.jpeg")

ggplot(UKB_CP_MDD_inflam_covariates_temp, aes(x=as.factor(recurrent_depression), y= CRP_log, group = recurrent_depression)) + 
  geom_boxplot(fill = c("#ffffff", "#FFE503"),alpha = 0.5) +
  theme_classic(base_size = 20)  +
  labs(title = "", x = "Probable Recurrent MDD Status", y = "Serum CRP (log, mg/l)")
dev.off()

## `Plot CRP in chronic pain depression comorbidity
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status <- factor(UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status  , levels=c("No Probable Recurrent MDD + No Chronic Pain",  "Probable Recurrent MDD + No Chronic Pain", "No Probable Recurrent MDD + Chronic Pain",
                                                                                                                                    "Probable Recurrent MDD + Chronic Pain"))
## Create temp df with no NAs in chronic pain status and CRP columns
UKB_CP_MDD_inflam_covariates_temp <- UKB_CP_MDD_inflam_covariates[!is.na(UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status) & !is.na(UKB_CP_MDD_inflam_covariates$CRP_log), ]


jpeg("~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/phenotypic/plots/CRP_CPMDD.jpeg")


ggplot(data=UKB_CP_MDD_inflam_covariates_temp) +
  (aes(x=as.factor(CP_MDD_comorbidity_status), y=CRP_log)) +
  geom_boxplot(fill = c("#ffffff", "#FFE503", "#2986cc", "#2ABA00"), alpha = 0.5) +
  theme_classic(base_size = 15) +
  labs(title = "", x = "Comorbidity Group", y = "Serum CRP (log, mg/l)") +
  scale_x_discrete(labels=c("No Probable Recurrent MDD + No Chronic Pain" = "CP-MDD-",
                            "Probable Recurrent MDD + No Chronic Pain" = "CP-MDD+",
                            "No Probable Recurrent MDD + Chronic Pain" = "CP+MDD-",
                            "Probable Recurrent MDD + Chronic Pain" = "CP+MDD+"))

dev.off()


## Plot sex interaction on association between CRP and chronic pain 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create temp df with no NAs in chronic pain status, CRP and gender columns
UKB_CP_MDD_inflam_covariates_temp <- UKB_CP_MDD_inflam_covariates[!is.na(UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status) & !is.na(UKB_CP_MDD_inflam_covariates$CRP_log) & !is.na(UKB_CP_MDD_inflam_covariates$sex), ]

## Remove those without CRP

## Recode recurrent MDD status column for plotting
UKB_CP_MDD_inflam_covariates_temp$chronic_pain_status <- dplyr::recode(UKB_CP_MDD_inflam_covariates_temp$chronic_pain_status, 
                                                                       `1` = "Chronic Pain",
                                                                       `0` = "No Chronic Pain")


UKB_CP_MDD_inflam_covariates_temp$chronic_pain_status <- factor(UKB_CP_MDD_inflam_covariates_temp$chronic_pain_status, levels = c("No Chronic Pain", "Chronic Pain"))

jpeg('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/phenotypic/plots/CP_status_CRP_sex_interaction.jpg', width = 700, height = 600)

UKB_CP_MDD_inflam_covariates_temp %>% 
  ggplot() +
  aes(x = chronic_pain_status, color = sex, group = sex, y = CRP_log) +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line") +
  theme_classic(base_size = 20) +
  ylab("Serum CRP (log, mg/l)") + 
  xlab("Chronic Pain Status") +
  scale_color_manual(values=c("#BA0000", "#2b83ba")) +
  scale_x_discrete(labels=c("0" = "No Chronic Pain",
                            "1" = "Chronic Pain"))

dev.off()


## Plot sex interaction on association between CRP and chronic pain group
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Remove those with no chronic pain site info for plotting and store as temporary df
UKB_CP_MDD_inflam_covariates_temp <- UKB_CP_MDD_inflam_covariates[!is.na(UKB_CP_MDD_inflam_covariates$chronic_pain_group) & !is.na(UKB_CP_MDD_inflam_covariates$sex), ]



UKB_CP_MDD_inflam_covariates_temp$chronic_pain_group <- factor(UKB_CP_MDD_inflam_covariates_temp$chronic_pain_group, levels = c("No Sites",
                                                                                                                                "1-2 Sites",
                                                                                                                                "3-4 Sites",
                                                                                                                                "5-7 Sites",
                                                                                                                                "Widespread"))

jpeg('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/phenotypic/plots/CP_group_CRP_sex_interaction.jpg', width = 650, height = 400)

UKB_CP_MDD_inflam_covariates_temp%>% 
  ggplot() +
  aes(x = chronic_pain_group, color = sex, group = sex, y = CRP_log) +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line") +
  theme_classic(base_size = 20) + 
  ylab("Serum CRP (log, mg/l)") + 
  xlab("Chronic Pain Group") +
  scale_color_manual(values=c("#BA0000", "#2b83ba"))


dev.off()

## Plot sex interaction on association between CRP and MDD
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Remove those with no chronic pain site info for plotting and store as temporary df
UKB_CP_MDD_inflam_covariates_temp <- UKB_CP_MDD_inflam_covariates[!is.na(UKB_CP_MDD_inflam_covariates$recurrent_depression) & !is.na(UKB_CP_MDD_inflam_covariates$sex), ]

jpeg('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/phenotypic/plots/MDD_status_CRP_sex_interaction.jpg', width = 700, height = 600)

UKB_CP_MDD_inflam_covariates_temp %>% 
  ggplot() +
  aes(x = recurrent_depression, color = sex, group = sex, y = CRP_log) +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line") +
  theme_classic(base_size = 20) +
  ylab("Serum CRP (log, mg/l)") + 
  xlab("Recurrent MDD Status") +
  scale_color_manual(values=c("#BA0000", "#2b83ba"))

dev.off()

## Plot sex interaction on assocation of CP+MDD and CRP
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Remove those with no CPMDD status for plotting and store as temporary df
UKB_CP_MDD_inflam_covariates_temp <- UKB_CP_MDD_inflam_covariates[!is.na(UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status),]


UKB_CP_MDD_inflam_covariates_temp$CP_MDD_comorbidity_status <- factor(UKB_CP_MDD_inflam_covariates_temp$CP_MDD_comorbidity_status, levels = c("No Probable Recurrent MDD + No Chronic Pain",
                                                                                                                                              "Probable Recurrent MDD + No Chronic Pain",
                                                                                                                                              "No Probable Recurrent MDD + Chronic Pain",
                                                                                                                                              "Probable Recurrent MDD + Chronic Pain"))

jpeg('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/phenotypic/plots/CPMDD_group_CRP_sex_interaction.jpg', width = 500, height = 400)


UKB_CP_MDD_inflam_covariates_temp %>% 
  ggplot() +
  aes(x = CP_MDD_comorbidity_status, color = sex, group = sex, y = CRP_log) +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line") +
  theme_classic(base_size = 13) + 
  ylab("Serum CRP (log, mg/l)") + 
  xlab("Comorbid Chronci Pain and Depression Group")+ 
  scale_color_manual(values=c("#BA0000", "#2b83ba")) +
  scale_x_discrete(labels=c("No Probable Recurrent MDD + No Chronic Pain" = "CP-MDD-",
                            "Probable Recurrent MDD + No Chronic Pain" = "CP-MDD+",
                            "No Probable Recurrent MDD + Chronic Pain" = "CP+MDD-",
                            "Probable Recurrent MDD + Chronic Pain" = "CP+MDD+"))

dev.off()



## Statistical analysis: interaction analysis (CRP ~ Chronic Pain Status + sex + sex:Chronic Pain Status + Covariates)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CP_status_CRP_sex_interaction_all_covar <- lm(CRP_log ~ chronic_pain_status * sex + alcohol_consumption + smoking_status + age + BMI_cat + deprivation,
                                                data=UKB_CP_MDD_inflam_covariates)

CRP_results <- data_frame(analysis = rep(NA, 100), beta = rep(NA, 100), sd = rep(NA, 100), p = rep(NA, 100), p.adjust = rep(NA, 100), nobs = rep(NA, 100))
CRP_results[1,1] <- "CP_sex"
CRP_results[1,2:4] <- t(summary(CP_status_CRP_sex_interaction_all_covar)$coefficient["chronic_pain_status:sexMale", c("Estimate", "Std. Error", "Pr(>|t|)")])
CRP_results[1,6] <- nobs(CP_status_CRP_sex_interaction_all_covar)
## Apply multiple comparison testing, n = 2 (accounts for two inflammatory biomarkers being used)
CRP_results[1,5] <- p.adjust(CRP_results[1,4], n = 2, method = "bonferroni")

## Statistical analysis: interaction analysis (CRP ~ Chronic Pain Site Groups + sex + sex:Chronic Pain Site Groups  + Covariates)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Set "no pain sites" as reference category
UKB_CP_MDD_inflam_covariates$chronic_pain_group <-  relevel(as.factor(UKB_CP_MDD_inflam_covariates$chronic_pain_group), ref = "No Sites")

CP_group_CRP_sex_interaction_all_covar <- lm(CRP_log ~ chronic_pain_group * sex + alcohol_consumption + smoking_status + age  + BMI_cat + deprivation, 
                                               data=UKB_CP_MDD_inflam_covariates)

summary(CP_group_CRP_sex_interaction_all_covar)

CRP_results[2:5,1] <- (gsub("chronic_pain_group", "sex cp group ", row.names(summary(CP_group_CRP_sex_interaction_all_covar)$coefficient)[14:17]))
CRP_results[2:5,2:4] <- (summary(CP_group_CRP_sex_interaction_all_covar)$coefficient[14:17, c("Estimate", "Std. Error", "Pr(>|t|)")])
CRP_results[2:5,6] <- nobs(CP_group_CRP_sex_interaction_all_covar)

## Apply multiple comparison testing, n = 8 (accounts for two inflammatory biomarkers being used)
CRP_results[2:5,5] <- (p.adjust(tail(summary(CP_group_CRP_sex_interaction_all_covar)$coefficients, n = 4)[,"Pr(>|t|)"], method = "bonferroni", n = 8))

## Statistical analysis: interaction analysis (CRP ~ Recurrent MDD Status + sex + sex:Recurrent MDD Status  + Covariates)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MDD_status_CRP_sex_interaction_all_covar <- lm(CRP_log ~ recurrent_depression * sex + alcohol_consumption + smoking_status + age + BMI_cat + deprivation, 
                                                 data=UKB_CP_MDD_inflam_covariates)
CRP_results[6,1] <- "MDD_sex"
CRP_results[6,2:4] <- t(summary(MDD_status_CRP_sex_interaction_all_covar)$coefficient["recurrent_depression:sexMale", c("Estimate", "Std. Error", "Pr(>|t|)")])
CRP_results[6,6] <- nobs(MDD_status_CRP_sex_interaction_all_covar)
## Apply multiple comparison testing, n = 2 (accounts for two inflammatory biomarkers being used)
CRP_results[6,5] <- p.adjust(CRP_results[1,4], n = 2, method = "bonferroni")

## Statistical analysis: interaction analysis (CRP ~ Comorbidity Group + sex + sex:Comorbidity Group  + Covariates)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Set "CP-MDD-" as reference category
UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status <- factor(UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status, 
                                                                 levels = c("No Probable Recurrent MDD + No Chronic Pain", 
                                                                            "Probable Recurrent MDD + No Chronic Pain", 
                                                                            "No Probable Recurrent MDD + Chronic Pain", 
                                                                            "Probable Recurrent MDD + Chronic Pain"))

UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status <- relevel(UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status, ref = "No Probable Recurrent MDD + No Chronic Pain")

CPMDD_status_CRP_sex_interaction_all_covar <- lm(CRP_log ~ CP_MDD_comorbidity_status * sex + alcohol_consumption + smoking_status + age + BMI_cat + deprivation,
                                                   data=UKB_CP_MDD_inflam_covariates)

summary(CPMDD_status_CRP_sex_interaction_all_covar)

CRP_results[7:9,1] <- (gsub("CP_MDD_comorbidity_status", "sex CPMDD group ", row.names(summary(CPMDD_status_CRP_sex_interaction_all_covar)$coefficient)[13:15]))
CRP_results[7:9,2:4] <- (summary(CPMDD_status_CRP_sex_interaction_all_covar)$coefficient[13:15, c("Estimate", "Std. Error", "Pr(>|t|)")])
CRP_results[7:9,6] <- nobs(CPMDD_status_CRP_sex_interaction_all_covar)

## Apply multiple comparison testing, n = 8 (accounts for two inflammatory biomarkers being used)
CRP_results[7:9,5] <- (p.adjust(tail(summary(CPMDD_status_CRP_sex_interaction_all_covar)$coefficients, n = 3)[,"Pr(>|t|)"], method = "bonferroni", n = 6))


## Statistical analysis: (x ~ covariate + group) - CRP ~ Chronic Pain Status 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create temp df with no NAs in chronic pain status and CRP columns
UKB_CP_MDD_inflam_covariates_temp <- UKB_CP_MDD_inflam_covariates[!is.na(UKB_CP_MDD_inflam_covariates$chronic_pain_status) & !is.na(UKB_CP_MDD_inflam_covariates$CRP_log), ]

UKB_CP_MDD_inflam_covariates_temp$chronic_pain_status <- as.factor(UKB_CP_MDD_inflam_covariates_temp$chronic_pain_status)

## run second model: adjusted for all covariates
model_2 <- lm(CRP_log ~ chronic_pain_status + BMI_cat + sex + alcohol_consumption + smoking_status + age + deprivation, data = UKB_CP_MDD_inflam_covariates_temp)
summary(model_2)
nobs(model_2)

CRP_results[10,1] <- "CP"
CRP_results[10,2:4] <- t(summary(model_2)$coefficient["chronic_pain_status1", c("Estimate", "Std. Error", "Pr(>|t|)")])
CRP_results[10,6] <- nobs(model_2)
## Apply multiple comparison testing, n = 2 (accounts for two inflammatory biomarkers being used)
CRP_results[10,5] <- p.adjust(CRP_results[10,4], n = 2, method = "bonferroni")


## Statistical analysis: (x ~ covariate + group) - CRP ~ Chronic Pain Group
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Change ranked and ordinal coviariates in order to run tuckey's (only accepts continuous covariates as ranked factors)
UKB_CP_MDD_inflam_covariates$age_rank <- as.factor(UKB_CP_MDD_inflam_covariates$age)
UKB_CP_MDD_inflam_covariates$deprivation_rank <- as.factor(UKB_CP_MDD_inflam_covariates$deprivation)
UKB_CP_MDD_inflam_covariates$smoking_status <- as.factor(UKB_CP_MDD_inflam_covariates$smoking_status)
UKB_CP_MDD_inflam_covariates$alcohol_consumption <- as.factor(UKB_CP_MDD_inflam_covariates$alcohol_consumption)


## Create temp df with no NAs in chronic pain status and CRP columns
UKB_CP_MDD_inflam_covariates_temp <- UKB_CP_MDD_inflam_covariates[!is.na(UKB_CP_MDD_inflam_covariates$chronic_pain_group) & !is.na(UKB_CP_MDD_inflam_covariates$CRP_log), ]

## run second model: adjusted for all covariates
lm1 <- lm(CRP_log ~ chronic_pain_group + BMI_cat + sex + alcohol_consumption + smoking_status + age + deprivation, data = UKB_CP_MDD_inflam_covariates_temp)
anova1 <- aov(CRP_log ~ chronic_pain_group + BMI_cat + sex + alcohol_consumption + smoking_status + age + deprivation, data = UKB_CP_MDD_inflam_covariates_temp)
summary(anova1)
## Post-hoc analysis
pairwise1 <- emmeans(lm1, list(pairwise ~ chronic_pain_group), adjust = "none")

CRP_results[11:20,1] <- as.data.frame(pairwise1$`pairwise differences of chronic_pain_group`)[,1]
CRP_results[11:20,2:4] <- as.data.frame(pairwise1$`pairwise differences of chronic_pain_group`)[,c("estimate", "SE", "p.value")]
CRP_results[11:20,6] <- nobs(lm1)
## Apply multiple comparison testing, n = 20 (accounts for two inflammatory biomarkers being used)
CRP_results[11:20,5] <- p.adjust(CRP_results$p[11:20], n = 20, method = "bonferroni")

## Statistical analysis: (x ~ covariate + group) - CRP ~ Depression Status 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create temp df with no NAs in chronic pain status and CRP columns
UKB_CP_MDD_inflam_covariates_temp <- UKB_CP_MDD_inflam_covariates[!is.na(UKB_CP_MDD_inflam_covariates$recurrent_depression) & !is.na(UKB_CP_MDD_inflam_covariates$CRP_log), ]

## Run linear regression
UKB_CP_MDD_inflam_covariates_temp$recurrent_depression <- as.factor(UKB_CP_MDD_inflam_covariates_temp$recurrent_depression)

## run second model: adjusted for all covariates
model_2 <- lm(CRP_log ~ recurrent_depression + BMI_cat + sex + alcohol_consumption + smoking_status + age + deprivation, data = UKB_CP_MDD_inflam_covariates_temp)
summary(model_2)

CRP_results[21,1] <- "MDD"
CRP_results[21,2:4] <- t(summary(model_2)$coefficient["recurrent_depression1", c("Estimate", "Std. Error", "Pr(>|t|)")])
CRP_results[21,6] <- nobs(model_2)
## Apply multiple comparison testing, n = 2 (accounts for two inflammatory biomarkers being used)
CRP_results[21,5] <- p.adjust(CRP_results[21,4], n = 2, method = "bonferroni")


## Statistical analysis: (x ~ covariate + group) - CRP ~ Comorbidity Group
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create temp df with no NAs in comorbidity group and CRP columns
UKB_CP_MDD_inflam_covariates_temp <- UKB_CP_MDD_inflam_covariates[!is.na(UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status) & !is.na(UKB_CP_MDD_inflam_covariates$CRP_log), ]

## run second model: adjusted for all covariates
lm1 <- lm(CRP_log ~ CP_MDD_comorbidity_status + BMI_cat + sex + alcohol_consumption + smoking_status + age + deprivation, data = UKB_CP_MDD_inflam_covariates_temp)
anova1 <- aov(CRP_log ~ CP_MDD_comorbidity_status + BMI_cat + sex + alcohol_consumption + smoking_status + age + deprivation, data = UKB_CP_MDD_inflam_covariates_temp)
summary(anova1)
## Post-hoc analysis
pairwise1 <- emmeans(lm1, list(pairwise ~ CP_MDD_comorbidity_status), adjust = "none")

CRP_results[22:27,1] <- as.data.frame(pairwise1$`pairwise differences of CP_MDD_comorbidity_status`)[,1]
CRP_results[22:27,2:4] <- as.data.frame(pairwise1$`pairwise differences of CP_MDD_comorbidity_status`)[,c("estimate", "SE", "p.value")]
CRP_results[22:27,6] <- nobs(lm1)
## Apply multiple comparison testing, n = 12 (accounts for two inflammatory biomarkers being used)
CRP_results[22:27,5] <- p.adjust(CRP_results$p[22:27], n = 12, method = "bonferroni")

## Save results
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(CRP_results, "~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/phenotypic/CRP_results.csv", quote = FALSE, row.names = FALSE)


