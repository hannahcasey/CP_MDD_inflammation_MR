## Load libraries

library(ggplot2)
library(patchwork)


## Run CRP and glycA association scripts to save plots in R environment
jpeg('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/phenotypic/plots/all_associations.jpg', width = 1300, height = 1800)

all_association <- ((CP_status_CRP | CP_status_glycA) / (CP_groups_CRP | CP_group_glycA) / (MDD_CRP | MDD_glycA) / (CPMDD_CRP | CPMDD_glycA)) + 
  plot_annotation(tag_levels = 'A') &
  theme(axis.title = element_text(size = 18))
all_association

dev.off()



jpeg('~/Desktop/PhD/projects/CP_MDD_inflammation_MR/output/phenotypic/plots/all_sex_interaction.jpg', width = 1300, height = 1700)

all_sex_interaction <- (CP_status_CRP_sex_interaction | CP_status_glycA_sex_interaction) / (CP_group_CRP_sex_interaction | CP_group_glycA_sex_interaction) /(MDD_CRP_sex_interaction | MDD_glycA_sex_interaction) + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect") &
  theme(axis.title = element_text(size = 18))
all_sex_interaction

dev.off()
