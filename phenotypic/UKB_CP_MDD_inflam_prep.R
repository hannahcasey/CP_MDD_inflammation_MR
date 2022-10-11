## Load libraries
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(Hmisc)

## Load in data
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Chronic pain phenotype - KJA Johnston
## Depression phenotype - N Sangha
## Alcohol intake frequency (1558)
## Smoking status (1239)
## Body mass index (21001)
## Townsend deprivation index at recruitment (189)
## Sex (31)
## Age (21003)
## CRP levels (30710)
## gkycA levels ()

## Load in chronic pain phenotype
UKB_chronic_pain_phenotype <- read.csv("~/Desktop/PhD/projects/UKB_chronic_pain_phenotype/output/UKB_chronic_pain_phenotype.csv")

##Quantify chronic pain sites
UKB_chronic_pain_phenotype$chronic_pain_group<- NA
UKB_chronic_pain_phenotype$chronic_pain_group[UKB_chronic_pain_phenotype$chronic_pain_sites == 0] <- "No Sites"
UKB_chronic_pain_phenotype$chronic_pain_group[UKB_chronic_pain_phenotype$chronic_pain_sites > 0 & UKB_chronic_pain_phenotype$chronic_pain_sites < 3] <- "1-2 Sites"
UKB_chronic_pain_phenotype$chronic_pain_group[UKB_chronic_pain_phenotype$chronic_pain_sites > 2 & UKB_chronic_pain_phenotype$chronic_pain_sites < 5] <- "3-4 Sites"
UKB_chronic_pain_phenotype$chronic_pain_group[UKB_chronic_pain_phenotype$chronic_pain_sites > 4 & UKB_chronic_pain_phenotype$chronic_pain_sites < 8] <- "5-7 Sites"
UKB_chronic_pain_phenotype$chronic_pain_group[UKB_chronic_pain_phenotype$chronic_pain_sites == 8] <- "Widespread"

## keep only ID, chronic pain sites, chronic pain status data
UKB_chronic_pain_phenotype <- UKB_chronic_pain_phenotype %>%
  select(n_eid, chronic_pain_sites, chronic_pain_status, chronic_pain_group)

## Load in depression phenotype
UKB_depression_phenotype <- read.csv("~/Desktop/PhD/projects/UKB_depression_phenotype/output/UKB_depression_phenotype.csv", header = T)
## Rename ID column and keep only ID and patient group data
UKB_depression_phenotype <- UKB_depression_phenotype %>%
  rename(n_eid = f_eid) %>%
  select(n_eid, patient_group, recurrent_depression, single_depression, no_depression)
## Add extra column indicating depression status
UKB_depression_phenotype$depression_status <- NA
UKB_depression_phenotype$depression_status[UKB_depression_phenotype$patient_group == "recurrent depression"] <- "Probable Recurrent MDD"
UKB_depression_phenotype$depression_status[UKB_depression_phenotype$patient_group == "single depression"] <- "Probable Single MDD"
UKB_depression_phenotype$depression_status[UKB_depression_phenotype$patient_group %in% c("control", "bipolar depression", "unipolar mania")] <- "No Probable MDD"


## Load in touchscreen questionnaire - alcohol, smoking
UKB_touchscreen <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2021-04-phenotypes-ukb44797/Touchscreen.rds")
## Extract alcohol and  smoking data
UKB_alcohol_smoking <- cbind(UKB_touchscreen[,1], UKB_touchscreen[grepl("f.1558.", names(UKB_touchscreen))],
                                  UKB_touchscreen[grepl("f.1239.", names(UKB_touchscreen))])
## Rename and keep ID, alcohol and smoking data
UKB_alcohol_smoking <- UKB_alcohol_smoking %>%
  rename(n_eid = `UKB_touchscreen[, 1]`,
         alcohol_consumption = f.1558.0.0,
         smoking_status = f.1239.0.0) %>%
  select(n_eid, alcohol_consumption, smoking_status)
## Remove redundant dataframe
rm(UKB_touchscreen)


## Load in physical measures data - BMI
UKB_physical_measures <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2018-10-phenotypes-ukb24262/PhysicalMeasures.rds")
## Extract ID column and BMI measurement
UKB_BMI <- cbind(UKB_physical_measures[,1], UKB_physical_measures[grepl("f.21001.", names(UKB_physical_measures))])
## Rename BMI column
UKB_BMI <- UKB_BMI %>%
  rename(n_eid = `UKB_physical_measures[, 1]`,
          BMI = f.21001.0.0) %>%
  select(n_eid, BMI)
## Remove redundant dataframes 
rm(UKB_physical_measures)


## Load in baseline characteristics - deprivation, sex
UKB_baseline_characteristics <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2021-04-phenotypes-ukb44797/BaselineCharacteristics.rds")
## Extract deprivation and sex daat
UKB_sex_deprivation <- cbind(UKB_baseline_characteristics[,1],
                                               UKB_baseline_characteristics[grepl("f.31.0.0", names(UKB_baseline_characteristics))],
                                               UKB_baseline_characteristics[grepl("f.189.0.0", names(UKB_baseline_characteristics))])
## Renames ID, sex and deprication column
UKB_sex_deprivation <- UKB_sex_deprivation %>%
  rename(n_eid = `UKB_baseline_characteristics[, 1]`,
         sex = `f.31.0.0`,
         deprivation = `f.189.0.0`)
## Remove redundant dataframes 
rm(UKB_baseline_characteristics)


## Load in recruitment data - age
UKB_recruitment <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2021-04-phenotypes-ukb44797/Recruitment.rds")
UKB_age <- cbind(UKB_recruitment[,1],                               UKB_recruitment[grepl("f.21003.0.0", names(UKB_recruitment))])
## rename ID and age columns
UKB_age <- UKB_age %>%
  rename(n_eid = `UKB_recruitment[, 1]`,
         age = f.21003.0.0)
## Remove redundant dataframes 
rm(UKB_recruitment)


## Load in assay data
UKB_assay <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2020-09-phenotypes-ukb43743/AssayResults.rds")
## f.30710 = C-reactive protein
## Extract ID column and columns pertaining to CRP measurement
UKB_CRP <- cbind(UKB_assay[,1], UKB_assay[grepl("f.30710.0.0", names(UKB_assay))])
## rename ID and age columns
UKB_CRP <- UKB_CRP %>%
  rename(n_eid = `UKB_assay[, 1]`,
         CRP = f.30710.0.0) %>%
  mutate(CRP_log = log(CRP))
## Remove redundant dataframes 
rm(UKB_assay)

## Load in NMR metabolomics data
UKB_NMR <- read.table("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2021-10-nmr-metabolomics-ukb48936/NMRMetabolomics.tsv.gz", header = TRUE)
## f.30710 = C-reactive protein
## Extract ID column and columns pretaining to glycA measurment
UKB_NMR_glycA <- cbind(UKB_NMR[,1], UKB_NMR[grepl("f.23480.", names(UKB_NMR))])

## Remove redundant dataframes 
rm(UKB_NMR)
## Rename ID column
UKB_NMR_glycA <- UKB_NMR_glycA %>%
  rename(n_eid = "UKB_NMR[, 1]",
         glycA_1 = f.23480.0.0,
         glycA_2 = f.23480.1.0)


## combine all filtered dataframes
df <- dplyr::left_join(UKB_chronic_pain_phenotype, UKB_depression_phenotype, by = "n_eid")
df <- left_join(df, UKB_sex_deprivation, by = "n_eid")
df <- left_join(df, UKB_alcohol_smoking, by = "n_eid")
df <- left_join(df, UKB_BMI, by = "n_eid")
df <- left_join(df, UKB_age, by = "n_eid")
df <- left_join(df, UKB_CRP, by = "n_eid")
df <- left_join(df, UKB_NMR_glycA, by = "n_eid")

UKB_CP_MDD_inflam_covariates <- df

## Replace CRP levels greater than 10mg/L with NA in first instance
UKB_CP_MDD_inflam_covariates$CRP[UKB_CP_MDD_inflam_covariates$CRP > 10] <- NA

## Indicate CP+MDD+ comorbidity status for those with available chronic pain and depression data
UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status <- NA
UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status[UKB_CP_MDD_inflam_covariates$recurrent_depression == 0  & UKB_CP_MDD_inflam_covariates$chronic_pain_status == 0] <- "No Probable Recurrent MDD + No Chronic Pain"
UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status[UKB_CP_MDD_inflam_covariates$recurrent_depression == 1  & UKB_CP_MDD_inflam_covariates$chronic_pain_status == 0] <- "Probable Recurrent MDD + No Chronic Pain"
UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status[UKB_CP_MDD_inflam_covariates$recurrent_depression == 0  & UKB_CP_MDD_inflam_covariates$chronic_pain_status == 1] <- "No Probable Recurrent MDD + Chronic Pain"
UKB_CP_MDD_inflam_covariates$CP_MDD_comorbidity_status[UKB_CP_MDD_inflam_covariates$recurrent_depression == 1 & UKB_CP_MDD_inflam_covariates$chronic_pain_status == 1] <- "Probable Recurrent MDD + Chronic Pain"

## Reformat data and convert continuous variables into ranked factors to run Tukey's post-hoc (continous variables not compatible)
## Reformat alcohol consumption data

## Daily or almost daily = 5
## Three or four times a week = 4
## Once or twice a week = 3
## One to three times a month = 2
## Special occasions only = 1
## Never = 0
## Prefer not to answer - NA

UKB_CP_MDD_inflam_covariates$alcohol_consumption <- as.character(UKB_CP_MDD_inflam_covariates$alcohol_consumption )

UKB_CP_MDD_inflam_covariates$alcohol_consumption[UKB_CP_MDD_inflam_covariates$alcohol_consumption == "Daily or almost daily"] <- 5
UKB_CP_MDD_inflam_covariates$alcohol_consumption[UKB_CP_MDD_inflam_covariates$alcohol_consumption == "Three or four times a week"] <- 4
UKB_CP_MDD_inflam_covariates$alcohol_consumption[UKB_CP_MDD_inflam_covariates$alcohol_consumption == "Once or twice a week"] <- 3
UKB_CP_MDD_inflam_covariates$alcohol_consumption[UKB_CP_MDD_inflam_covariates$alcohol_consumption == "One to three times a month"] <- 2
UKB_CP_MDD_inflam_covariates$alcohol_consumption[UKB_CP_MDD_inflam_covariates$alcohol_consumption == "Special occasions only"] <- 1
UKB_CP_MDD_inflam_covariates$alcohol_consumption[UKB_CP_MDD_inflam_covariates$alcohol_consumption == "Never"] <- 0
UKB_CP_MDD_inflam_covariates$alcohol_consumption[UKB_CP_MDD_inflam_covariates$alcohol_consumption == "Prefer not to answer"] <- NA

unique(UKB_CP_MDD_inflam_covariates$alcohol_consumption)
UKB_CP_MDD_inflam_covariates$alcohol_consumption <- as.factor(UKB_CP_MDD_inflam_covariates$alcohol_consumption)                                  

## Reformat smoking data
## Yes, on most or all days = 2
## Only occasionally = 1
## No = 0
## Prefer not to answer - NA

UKB_CP_MDD_inflam_covariates$smoking_status <- as.character(UKB_CP_MDD_inflam_covariates$smoking_status )

UKB_CP_MDD_inflam_covariates$smoking_status[UKB_CP_MDD_inflam_covariates$smoking_status == "Yes, on most or all days"] <- 2
UKB_CP_MDD_inflam_covariates$smoking_status[UKB_CP_MDD_inflam_covariates$smoking_status == "Only occasionally"] <- 1
UKB_CP_MDD_inflam_covariates$smoking_status[UKB_CP_MDD_inflam_covariates$smoking_status == "No"] <- 0
UKB_CP_MDD_inflam_covariates$smoking_status[UKB_CP_MDD_inflam_covariates$smoking_status == "Prefer not to answer"] <- NA

unique(UKB_CP_MDD_inflam_covariates$smoking_status)
UKB_CP_MDD_inflam_covariates$smoking_status <- as.factor(UKB_CP_MDD_inflam_covariates$smoking_status)

## Reformat BMI data
## <18.5 = 0
## 18.5 - 24.9 = 1
## 25 - 29.9 = 2
## >= 30 = 3

UKB_CP_MDD_inflam_covariates$BMI_cat <- NA
UKB_CP_MDD_inflam_covariates$BMI_cat[UKB_CP_MDD_inflam_covariates$BMI < 18.5] <- "Underweight"
UKB_CP_MDD_inflam_covariates$BMI_cat[(UKB_CP_MDD_inflam_covariates$BMI > 18.5 & UKB_CP_MDD_inflam_covariates$BMI < 25)] <- "Normal"
UKB_CP_MDD_inflam_covariates$BMI_cat[(UKB_CP_MDD_inflam_covariates$BMI >= 25 & UKB_CP_MDD_inflam_covariates$BMI < 30)] <- "Overweight"
UKB_CP_MDD_inflam_covariates$BMI_cat[UKB_CP_MDD_inflam_covariates$BMI >= 30] <- "Obese"

unique(UKB_CP_MDD_inflam_covariates$BMI_cat)
UKB_CP_MDD_inflam_covariates$BMI_cat <- as.factor(UKB_CP_MDD_inflam_covariates$BMI_cat)


## Rank Deprivation score
UKB_CP_MDD_inflam_covariates$deprivation_rank <- as.numeric(cut2(UKB_CP_MDD_inflam_covariates$deprivation, g = 10))
UKB_CP_MDD_inflam_covariates$deprivation_rank <- as.factor(UKB_CP_MDD_inflam_covariates$deprivation_rank)


## Rank age data
UKB_CP_MDD_inflam_covariates$age_rank <- as.numeric(cut2(UKB_CP_MDD_inflam_covariates$age, g = 10))
UKB_CP_MDD_inflam_covariates$age_rank <- as.factor(UKB_CP_MDD_inflam_covariates$age_rank)


## Rank CRP data
## <1 = 0
## 1=3 = 1
## 3.01 - 10 = 2
## >10 = 3
UKB_CP_MDD_inflam_covariates$CRP_cat <- NA
UKB_CP_MDD_inflam_covariates$CRP_cat[UKB_CP_MDD_inflam_covariates$CRP < 1] <- 0
UKB_CP_MDD_inflam_covariates$CRP_cat[(UKB_CP_MDD_inflam_covariates$CRP >= 1 & UKB_CP_MDD_inflam_covariates$CRP <= 3)] <- 1
UKB_CP_MDD_inflam_covariates$CRP_cat[(UKB_CP_MDD_inflam_covariates$CRP > 3 & UKB_CP_MDD_inflam_covariates$CRP <= 10)] <- 2
UKB_CP_MDD_inflam_covariates$CRP_cat[UKB_CP_MDD_inflam_covariates$CRP >= 10] <- 3

## Mutate principal components into ranked factors
UKB_CP_MDD_inflam_covariates <- UKB_CP_MDD_inflam_covariates %>%
  mutate_at(vars(contains('PC')), funs(rank = (as.numeric(cut2(., g = 10))))) %>%
  mutate_at(vars(contains('PC_rank')), funs(as.factor))

## Save formatted covariate data to resources subdirectory
write.csv(UKB_CP_MDD_inflam_covariates, "~/Desktop/PhD/projects/CP_MDD_inflammation_MR//resources/UKB_CP_MDD_inflam_covariates.csv", row.names = F)


