library(dplyr)

## Descriptive Statistics of CRP in UKB
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load in assay data
UKBAssayResults <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2020-09-phenotypes-ukb43743/AssayResults.rds")

## Limit to only IID and first instance CRP
UKBAssayResults_CRP <- UKBAssayResults[,c("f.eid", "f.30710.0.0")]
## Remove particpants with no recorded CRP (first instance)
UKBAssayResults_CRP <- UKBAssayResults_CRP[!is.na(UKBAssayResults_CRP$f.30710.0.0), ]
## Remove those with CRP greater than 10
UKBAssayResults_CRP_below10 <- UKBAssayResults_CRP[UKBAssayResults_CRP$f.30710.0.0 < 10,]

## Calculate number of participants in each CRP catagory 
length(na.omit(UKBAssayResults_CRP_below10$f.30710.0[UKBAssayResults_CRP_below10$f.30710.0 <= 1]))/length(na.omit(UKBAssayResults_CRP_below10$f.30710.0))
## 41% have less than 1mg/L
length(na.omit(UKBAssayResults_CRP_below10$f.30710.0[(UKBAssayResults_CRP_below10$f.30710.0 > 1 & UKBAssayResults_CRP_below10$f.30710.0 <= 3)]))/length(na.omit(UKBAssayResults_CRP_below10$f.30710.0))
## 39% have between 1 and 3 mg/L
length(na.omit(UKBAssayResults_CRP_below10$f.30710.0[(UKBAssayResults_CRP_below10$f.30710.0 > 3 & UKBAssayResults_CRP_below10$f.30710.0 <= 10)]))/length(na.omit(UKBAssayResults_CRP_below10$f.30710.0))
## 19 have between 3 and 10 mg/L

## Load in physical measures data and keep BMI
UKBPhysicalMeasures <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2018-10-phenotypes-ukb24262/PhysicalMeasures.rds")
UKBPhysicalMeasuresBMI <- cbind(UKBPhysicalMeasures[,1], UKBPhysicalMeasures[grepl("f.21001.", names(UKBPhysicalMeasures))])
## Rename ID column
names(UKBPhysicalMeasuresBMI)[names(UKBPhysicalMeasuresBMI) == "UKBPhysicalMeasures[, 1]"] <- "f.eid"

## Load in baseline characteristics and keep age, sex, deprivation
UKBBaselinecharacteristics <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2018-10-phenotypes-ukb24262/BaselineCharacteristics.rds")
UKBBaselinecharacteristicsAgeSexDep <- cbind(UKBBaselinecharacteristics[,1],
                                             UKBBaselinecharacteristics[grepl("f.31.0.0", names(UKBBaselinecharacteristics))],
                                             UKBBaselinecharacteristics[grepl("f.189.0.0", names(UKBBaselinecharacteristics))])
## Rename ID column
names(UKBBaselinecharacteristicsAgeSexDep)[names(UKBBaselinecharacteristicsAgeSexDep) == "UKBBaselinecharacteristics[, 1]"] <- "f.eid"


## Load in touchscreen questionnaire - alcohol, smoking
UKB_touchscreen <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2021-04-phenotypes-ukb44797/Touchscreen.rds")

## Extract alcohol and  smoking data
UKBAlcoholSmoking <- cbind(UKB_touchscreen[,1], UKB_touchscreen[grepl("f.1558.", names(UKB_touchscreen))],
                                  UKB_touchscreen[grepl("f.1239.", names(UKB_touchscreen))])

## Rename ID column
names(UKBAlcoholSmoking)[names(UKBAlcoholSmoking) == "UKB_touchscreen[, 1]"] <- "f.eid"

## Load in baseline characteristics and keep age, sex, deprivation
UKBRecruitment <- readRDS("/Volumes/GenScotDepression/data/ukb/phenotypes/fields/2018-10-phenotypes-ukb24262/Recruitment.rds")
UKBRecruitmentAge <- cbind(UKBRecruitment[,1],
                           UKBRecruitment[grepl("f.21003.0.0", names(UKBRecruitment))])
## Rename ID column
names(UKBRecruitmentAge)[names(UKBRecruitmentAge) == "UKBRecruitment[, 1]"] <- "f.eid"


## Join baseline, physical characteristics and touchscreen data to CRP data
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
UKB_CRP_BMI <- left_join(UKBAssayResults_CRP_below10, UKBPhysicalMeasuresBMI, by = "f.eid")
UKB_CRP_BMI_baseline <- left_join(UKB_CRP_BMI, UKBBaselinecharacteristicsAgeSexDep, by = "f.eid")
UKB_CRP_BMI_baseline_age <- left_join(UKB_CRP_BMI_baseline, UKBRecruitmentAge, by = "f.eid")
UKB_CRP_BMI_baseline_age_alcohol_smoking <- left_join(UKB_CRP_BMI_baseline_age, UKBAlcoholSmoking, by = "f.eid")


## Split up into CRP groups
All_below1 <- UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.30710.0.0 <= 1,]
All_between13 <- UKB_CRP_BMI_baseline_age_alcohol_smoking[(UKB_CRP_BMI_baseline_age_alcohol_smoking$f.30710.0.0 > 1 & UKB_CRP_BMI_baseline_age_alcohol_smoking$f.30710.0 <= 3),]
All_between310 <- UKB_CRP_BMI_baseline_age_alcohol_smoking[(UKB_CRP_BMI_baseline_age_alcohol_smoking$f.30710.0.0 > 3 & UKB_CRP_BMI_baseline_age_alcohol_smoking$f.30710.0 <= 10),]
## Make sure number of rows adds up to 449107
nrow(All_below1) + nrow(All_between13) + nrow(All_between310)


## Proportion of females in each CRP group
nrow(All_below1[All_below1$f.31.0.0 == "Female",]) / nrow(All_below1)
nrow(All_between13[All_between13$f.31.0.0 == "Female",]) / nrow(All_between13)
nrow(All_between310[All_between310$f.31.0.0 == "Female",]) / nrow(All_between310)

## Statistical analysis - do proportions of women in each group differ significantly
total_CRP_group <- c(nrow(All_below1), nrow(All_between13), nrow(All_between310))

res <- chisq.test(total_CRP_group, p = c(nrow(All_below1[All_below1$f.31.0.0 == "Female",]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.31.0.0 == "Female",]),
                                               nrow(All_between13[All_between13$f.31.0.0 == "Female",]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.31.0.0 == "Female",]),
                                               nrow(All_between310[All_between310$f.31.0.0 == "Female",]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.31.0.0 == "Female",])))
res$method


## Proportion of participants overweight/obese in each CRP group
nrow(All_below1[All_below1$f.21001.0.0 > 25,]) / nrow(All_below1)
nrow(All_between13[All_between13$f.21001.0.0 > 25,]) / nrow(All_between13)
nrow(All_between310[All_between310$f.21001.0.0 > 25,]) / nrow(All_between310)

## Statistical analysis - do proportions of obesity in each group differ significantly
res <- chisq.test(total_CRP_group, p = c(nrow(All_below1[All_below1$f.21001.0.0 > 25,]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.21001.0.0 > 25,]),
                                         nrow(All_between13[All_between13$f.21001.0.0 > 25,]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.21001.0.0 > 25,]),
                                         nrow(All_between310[All_between310$f.21001.0.0 > 25,]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.21001.0.0 > 25,])))
res

## Average age in each CRP group
mean(All_below1$f.21003.0.0)
mean(All_between13$f.21003.0.0)
mean(All_between310$f.21003.0.0)

## Statistical analysis - does age differ significantly between CRP groups?
UKB_CRP_BMI_baseline_age_alcohol_smoking$CRP_group <- NA
UKB_CRP_BMI_baseline_age_alcohol_smoking$CRP_group[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.30710.0.0 <= 1] <- "below_1"
UKB_CRP_BMI_baseline_age_alcohol_smoking$CRP_group[(UKB_CRP_BMI_baseline_age_alcohol_smoking$f.30710.0.0 > 1 & UKB_CRP_BMI_baseline_age_alcohol_smoking$f.30710.0 <= 3)] <- "between_1_3"
UKB_CRP_BMI_baseline_age_alcohol_smoking$CRP_group[(UKB_CRP_BMI_baseline_age_alcohol_smoking$f.30710.0.0 > 3 & UKB_CRP_BMI_baseline_age_alcohol_smoking$f.30710.0 <= 10)] <- "between_3_10"
freq_table(UKB_CRP_BMI_baseline_age_alcohol_smoking$CRP_group)

kruskal.test(f.21003.0.0 ~ CRP_group, data = UKB_CRP_BMI_baseline_age_alcohol_smoking)

## Average deprivation level in each CRP group
mean(na.omit(All_below1$f.189.0.0))
mean(na.omit(All_between13$f.189.0.0))
mean(na.omit(All_between310$f.189.0.0))

## Statistical analysis - do deprivation differ between CRP groups
hist(UKB_CRP_BMI_baseline_age_alcohol_smoking$f.189.0.0)
kruskal.test(f.189.0.0 ~ CRP_group, data = UKB_CRP_BMI_baseline_age_alcohol_smoking)

## Proportion of participants who consumes alcohol daily/ almost daily
nrow(All_below1[All_below1$f.1558.0.0 == "Female",]) / nrow(All_below1)
nrow(All_between13[All_between13$f.31.0.0 == "Female",]) / nrow(All_between13)
nrow(All_between310[All_between310$f.31.0.0 == "Female",]) / nrow(All_between310)

## Statistical analysis - do proportions of women in each group differ significantly
total_CRP_group <- c(nrow(All_below1), nrow(All_between13), nrow(All_between310))

res <- chisq.test(total_CRP_group, p = c(nrow(All_below1[All_below1$f.31.0.0 == "Female",]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.31.0.0 == "Female",]),
                                         nrow(All_between13[All_between13$f.31.0.0 == "Female",]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.31.0.0 == "Female",]),
                                         nrow(All_between310[All_between310$f.31.0.0 == "Female",]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.31.0.0 == "Female",])))
res




## Load in chronic pain and depresison datasets that have been merged with CRP data
UKB_CP_CRP_PM <- read.table("~/Desktop/PhD/projects/UKBChronicPainDepressionCRP/resources/UKB_CP_CRP_PM.txt", header = T)
UKB_MDD_CRP_PM <- read.table("~/Desktop/PhD/projects/UKBChronicPainDepressionCRP/resources/UKB_MDD_CRP_PM.txt", header = T, row.names = NULL)


## Split up into CRP groups
UKB_CP_CRP_PMbelow1 <- UKB_CP_CRP_PM[UKB_CP_CRP_PM$f.30710.0.0 <= 1,]
UKB_CP_CRP_PMbetween13 <- UKB_CP_CRP_PM[(UKB_CP_CRP_PM$f.30710.0.0 > 1 & UKB_CP_CRP_PM$f.30710.0 <= 3),]
UKB_CP_CRP_PMbetween310 <- UKB_CP_CRP_PM[(UKB_CP_CRP_PM$f.30710.0.0 > 3 & UKB_CP_CRP_PM$f.30710.0 <= 10),]

UKB_MDD_CRP_PMbelow1 <- UKB_MDD_CRP_PM[UKB_MDD_CRP_PM$f.30710.0.0 <= 1,]
UKB_MDD_CRP_PMbetween13 <- UKB_MDD_CRP_PM[(UKB_MDD_CRP_PM$f.30710.0.0 > 1 & UKB_MDD_CRP_PM$f.30710.0 <= 3),]
UKB_MDD_CRP_PMbetween310 <- UKB_MDD_CRP_PM[(UKB_MDD_CRP_PM$f.30710.0.0 > 3 & UKB_MDD_CRP_PM$f.30710.0 <= 10),]

## Ensure rows sum to 440705 and 150318
nrow(UKB_CP_CRP_PMbelow1) + nrow(UKB_CP_CRP_PMbetween13) + nrow(UKB_CP_CRP_PMbetween310)
nrow(UKB_MDD_CRP_PMbelow1) + nrow(UKB_MDD_CRP_PMbetween13) + nrow(UKB_MDD_CRP_PMbetween310)

##Proportion of chronic pain in each CRP group
nrow(UKB_CP_CRP_PMbelow1[UKB_CP_CRP_PMbelow1$ChronicPainStatus == 1,]) / nrow(UKB_CP_CRP_PMbelow1)
nrow(UKB_CP_CRP_PMbetween13[UKB_CP_CRP_PMbetween13$ChronicPainStatus == 1,]) / nrow(UKB_CP_CRP_PMbetween13)
nrow(UKB_CP_CRP_PMbetween310[UKB_CP_CRP_PMbetween310$ChronicPainStatus == 1,]) / nrow(UKB_CP_CRP_PMbetween310)


## Statistical analysis - do proportions of obesity in each group differ significantly
total_CRP_group <- c(nrow(UKB_CP_CRP_PMbelow1), nrow(UKB_CP_CRP_PMbetween13), nrow(UKB_CP_CRP_PMbetween310))

res <- chisq.test(total_CRP_group, p = c(nrow(UKB_CP_CRP_PMbelow1[UKB_CP_CRP_PMbelow1$ChronicPainStatus == 1,]) / nrow(UKB_CP_CRP_PM[UKB_CP_CRP_PM$ChronicPainStatus == 1,]),
                                         nrow(UKB_CP_CRP_PMbetween13[UKB_CP_CRP_PMbetween13$ChronicPainStatus == 1,]) / nrow(UKB_CP_CRP_PM[UKB_CP_CRP_PM$ChronicPainStatus == 1,]),
                                         nrow(UKB_CP_CRP_PMbetween310[UKB_CP_CRP_PMbetween310$ChronicPainStatus == 1,]) / nrow(UKB_CP_CRP_PM[UKB_CP_CRP_PM$ChronicPainStatus == 1,])))
res

##Proportion of depression in each CRP group
nrow(UKB_MDD_CRP_PMbelow1[UKB_MDD_CRP_PMbelow1$DepressionStatus == "Probable Recurrent MDD" ,]) / nrow(UKB_MDD_CRP_PMbelow1)
nrow(UKB_MDD_CRP_PMbetween13[UKB_MDD_CRP_PMbetween13$DepressionStatus == "Probable Recurrent MDD",]) / nrow(UKB_MDD_CRP_PMbetween13)
nrow(UKB_MDD_CRP_PMbetween310[UKB_MDD_CRP_PMbetween310$DepressionStatus == "Probable Recurrent MDD",]) / nrow(UKB_MDD_CRP_PMbetween310)


## Proportion of people who smoke on most or all days in each CRP group
nrow(All_below1[All_below1$f.1558.0.0 == "Daily or almost daily",]) / nrow(All_below1)
nrow(All_between13[All_between13$f.1558.0.0 == "Daily or almost daily",]) / nrow(All_between13)
nrow(All_between310[All_between310$f.1558.0.0 == "Daily or almost daily",]) / nrow(All_between310)


## Statistical analysis - do proportions of people who smoke on most or all days in each CRP group differ?
total_CRP_group <- c(nrow(All_below1), nrow(All_between13), nrow(All_between310))

res <- chisq.test(total_CRP_group, p = c(nrow(All_below1[All_below1$f.1558.0.0 == "Daily or almost daily",]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.1558.0.0 == "Daily or almost daily",]),
                                         nrow(All_between13[All_between13$f.1558.0.0 == "Daily or almost daily",]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.1558.0.0 == "Daily or almost daily",]),
                                         nrow(All_between310[All_between310$f.1558.0.0 == "Daily or almost daily",]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.1558.0.0 == "Daily or almost daily",])))
res

## Proportion of people who  who consume alcohol daily/ almost daily in each CRP group
nrow(All_below1[All_below1$f.1239.0.0 == "Yes, on most or all days",]) / nrow(All_below1)
nrow(All_between13[All_between13$f.1239.0.0 == "Yes, on most or all days",]) / nrow(All_between13)
nrow(All_between310[All_between310$f.1239.0.0 == "Yes, on most or all days",]) / nrow(All_between310)


## Statistical analysis - do proportions of people who consume alcohol daily/ almost daily in each CRP group differ?
total_CRP_group <- c(nrow(All_below1), nrow(All_between13), nrow(All_between310))

res <- chisq.test(total_CRP_group, p = c(nrow(All_below1[All_below1$f.1239.0.0 == "Yes, on most or all days",]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.1239.0.0 == "Yes, on most or all days",]),
                                         nrow(All_between13[All_between13$f.1239.0.0 == "Yes, on most or all days",]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.1239.0.0 == "Yes, on most or all days",]),
                                         nrow(All_between310[All_between310$f.1239.0.0 == "Yes, on most or all days",]) / nrow(UKB_CRP_BMI_baseline_age_alcohol_smoking[UKB_CRP_BMI_baseline_age_alcohol_smoking$f.1239.0.0 == "Yes, on most or all days",])))
res


## Change class of stratifiled CRP dataframes to dataframes

## Proportion of comorbid chronic pain and depression in each CRP group
## First create dataframe with both chronic pain and depression data
UKB_comorbid_CRP_PMbelow1 <- inner_join(UKB_MDD_CRP_PMbelow1, UKB_CP_CRP_PMbelow1, by = "n_eid")
UKB_comorbid_CRP_PMbetween13<- inner_join(UKB_MDD_CRP_PMbetween13, UKB_CP_CRP_PMbetween13, by = "n_eid")
UKB_comorbid_CRP_PMbetween310 <- inner_join(UKB_MDD_CRP_PMbetween310, UKB_CP_CRP_PMbetween310, by = "n_eid")

## Create new column to store comorbidity status
UKB_comorbid_CRP_PMbelow1$comorbid_status <- NA
UKB_comorbid_CRP_PMbelow1$comorbid_status[(UKB_comorbid_CRP_PMbelow1$RecurrentMDDStatus == "Probable Recurrent MDD" & UKB_comorbid_CRP_PMbelow1$ChronicPainStatus == 1)] <- 1
UKB_comorbid_CRP_PMbelow1$comorbid_status[(UKB_comorbid_CRP_PMbelow1$RecurrentMDDStatus == "Probable Recurrent MDD" & UKB_comorbid_CRP_PMbelow1$ChronicPainStatus == 0)] <- 0
UKB_comorbid_CRP_PMbelow1$comorbid_status[(UKB_comorbid_CRP_PMbelow1$RecurrentMDDStatus == "No Probable Recurrent MDD" & UKB_comorbid_CRP_PMbelow1$ChronicPainStatus == 1)] <- 0
UKB_comorbid_CRP_PMbelow1$comorbid_status[(UKB_comorbid_CRP_PMbelow1$RecurrentMDDStatus == "No Probable Recurrent MDD" & UKB_comorbid_CRP_PMbelow1$ChronicPainStatus == 0)] <- 0

UKB_comorbid_CRP_PMbetween13$comorbid_status <- NA
UKB_comorbid_CRP_PMbetween13$comorbid_status[(UKB_comorbid_CRP_PMbetween13$RecurrentMDDStatus == "Probable Recurrent MDD" & UKB_comorbid_CRP_PMbetween13$ChronicPainStatus == 1)] <- 1
UKB_comorbid_CRP_PMbetween13$comorbid_status[(UKB_comorbid_CRP_PMbetween13$RecurrentMDDStatus == "Probable Recurrent MDD" & UKB_comorbid_CRP_PMbetween13$ChronicPainStatus == 0)] <- 0
UKB_comorbid_CRP_PMbetween13$comorbid_status[(UKB_comorbid_CRP_PMbetween13$RecurrentMDDStatus == "No Probable Recurrent MDD" & UKB_comorbid_CRP_PMbetween13$ChronicPainStatus == 1)] <- 0
UKB_comorbid_CRP_PMbetween13$comorbid_status[(UKB_comorbid_CRP_PMbetween13$RecurrentMDDStatus == "No Probable Recurrent MDD" & UKB_comorbid_CRP_PMbetween13$ChronicPainStatus == 0)] <- 0

UKB_comorbid_CRP_PMbetween310$comorbid_status <- NA
UKB_comorbid_CRP_PMbetween310$comorbid_status[(UKB_comorbid_CRP_PMbetween310$RecurrentMDDStatus == "Probable Recurrent MDD" & UKB_comorbid_CRP_PMbetween310$ChronicPainStatus == 1)] <- 1
UKB_comorbid_CRP_PMbetween310$comorbid_status[(UKB_comorbid_CRP_PMbetween310$RecurrentMDDStatus == "Probable Recurrent MDD" & UKB_comorbid_CRP_PMbetween310$ChronicPainStatus == 0)] <- 0
UKB_comorbid_CRP_PMbetween310$comorbid_status[(UKB_comorbid_CRP_PMbetween310$RecurrentMDDStatus == "No Probable Recurrent MDD" & UKB_comorbid_CRP_PMbetween310$ChronicPainStatus == 1)] <- 0
UKB_comorbid_CRP_PMbetween310$comorbid_status[(UKB_comorbid_CRP_PMbetween310$RecurrentMDDStatus == "No Probable Recurrent MDD" & UKB_comorbid_CRP_PMbetween310$ChronicPainStatus == 0)] <- 0

## Delete rows with NA in comorbid status column 
UKB_comorbid_CRP_PMbelow1 <- UKB_comorbid_CRP_PMbelow1[!is.na(UKB_comorbid_CRP_PMbelow1$comorbid_status), ]
UKB_comorbid_CRP_PMbetween13 <- UKB_comorbid_CRP_PMbetween13[!is.na(UKB_comorbid_CRP_PMbetween13$comorbid_status), ]
UKB_comorbid_CRP_PMbetween310 <- UKB_comorbid_CRP_PMbetween310[!is.na(UKB_comorbid_CRP_PMbelow1$comorbid_status), ]

nrow(UKB_comorbid_CRP_PMbelow1[(UKB_comorbid_CRP_PMbelow1$comorbid_status == 1),]) / nrow(UKB_comorbid_CRP_PMbelow1)
nrow(UKB_comorbid_CRP_PMbetween13[(UKB_comorbid_CRP_PMbetween13$comorbid_status == 1),]) / nrow(UKB_comorbid_CRP_PMbetween13)
nrow(UKB_comorbid_CRP_PMbetween310[(UKB_comorbid_CRP_PMbetween310$comorbid_status == 1),]) / nrow(UKB_comorbid_CRP_PMbetween310)


## Statistical analysis - do proportions of obesity in each group differ significantly
total_CRP_group <- c(nrow(UKB_comorbid_CRP_PMbelow1), nrow(UKB_comorbid_CRP_PMbetween13), nrow(UKB_comorbid_CRP_PMbetween310))

res <- chisq.test(total_CRP_group, p = c(nrow(UKB_MDD_CRP_PMbelow1[UKB_MDD_CRP_PMbelow1$DepressionStatus == "Probable Recurrent MDD",]) / nrow(UKB_MDD_CRP_PM[UKB_MDD_CRP_PM$DepressionStatus == "Probable Recurrent MDD",]),
                                         nrow(UKB_MDD_CRP_PMbetween13[UKB_MDD_CRP_PMbetween13$DepressionStatus == "Probable Recurrent MDD",]) / nrow(UKB_MDD_CRP_PM[UKB_MDD_CRP_PM$DepressionStatus == "Probable Recurrent MDD",]),
                                         nrow(UKB_MDD_CRP_PMbetween310[UKB_MDD_CRP_PMbetween310$DepressionStatus == "Probable Recurrent MDD",]) / nrow(UKB_MDD_CRP_PM[UKB_MDD_CRP_PM$DepressionStatus == "Probable Recurrent MDD",])))
res
