##### Script F #####
##### R script to reproduce Figures Fig 8b and 8g in GO-2-Substrates paper #####
##### Note - run Script A1 first


### Generate ROC and PR curves and plot for Fig 8b and 8g
cutsites_subset_unique_fimo_rank = cutsites[match(unique(cutsites$Gene_name), cutsites$Gene_name),]

library(precrec)
#calculate ROC and precision-recall for sequence module
subset_fimo = cutsites_subset_unique_fimo_rank %>% filter(!is.na(cutsites_subset_unique_fimo_rank$outcome_test_data))
subset_fimo = subset(cutsites_subset_unique_fimo_rank, select = c(outcome_test_data, fimo_rank_normalised))

#note - have to reverse normalised rank here (1= best) precrec doesnt allow direction setting
subset_fimo$fimo_rank_normalised = 1-subset_fimo$fimo_rank_normalised

tested_indices_fimo = which(!is.na(subset_fimo$outcome_test_data))
msmdat_fimo <- mmdata(subset_fimo$fimo_rank_normalised[tested_indices_fimo], 
                      subset_fimo$outcome_test_data[tested_indices_fimo])
sscurves <- evalmod(msmdat_fimo)
# plot Fig 8b and 8g
pdf("Fig_8b_8g.pdf") 
# 2. Create a plot
plot(sscurves)
# Close the pdf file
dev.off()
