##### Script E #####
##### R script to reproduce Figures Fig 3a and b, Fig 4a, b and c, Fig 8c-e and h-i, 
##### Supplementary Fig 1, Supplementary Fig 2 in GO-2-Substrates paper
##### Note, run Script A first

## Append experimental outcomes onto cutsites
assay_candidates = read.csv("cotransfection_assay_outcomes.csv", header=TRUE)
cutsites$literature_cleaved = NA
cutsites$cleaved = NA
cutsites$cleaved_in_assay = NA
cutsites$tested_in_assay = NA
cutsites$tested_so_far = NA
cutsites$outcome = NA
cutsites$ever_tested_in_assay = NA
cutsites$outcome_test_data = NA
cutsites$sequence_cut = NA

for(i in 1:nrow(assay_candidates)) {
  for(j in 1:nrow(cutsites)) {
    if(identical(as.vector(assay_candidates$Gene_name[i]), as.vector(cutsites$Gene_name[j]))) {
      cutsites$ever_tested_in_assay[j] = assay_candidates$ever_tested_in_assay[i]
      cutsites$literature_cleaved[j] = assay_candidates$literature_cleaved[i] 
      cutsites$cleaved[j] = assay_candidates$cleaved[i]
      cutsites$cleaved_in_assay[j] = assay_candidates$cleaved_in_assay[i] 
      cutsites$tested_in_assay[j] = assay_candidates$tested_in_assay[i]
      cutsites$tested_so_far[j] = assay_candidates$tested_so_far[i]
      cutsites$outcome[j] = as.numeric(assay_candidates$outcome[i])
      cutsites$outcome_test_data[j] = as.numeric(assay_candidates$outcome_test_data[i])
      print(assay_candidates[i,1])
    }
  }
}

### Prepare data for visualisation

cutsites_subset_unique_GO2Substrates_rank = cutsites
cutsites_subset_unique_GO2Substrates_rank$gene_rank = rank(cutsites_subset_unique_GO2Substrates_rank$final_rank_normalised, ties.method = "min")
cutsites_subset_unique_GO2Substrates_rank$go2substrates_rank = cutsites_subset_unique_GO2Substrates_rank$gene_rank

# extract the first instance of each unique gene from cutsites
# based on sequence rank
cutsites = cutsites[order(cutsites$sequence_rank),]
cutsites_subset_unique_sequence_rank = cutsites[match(unique(cutsites$Entry.name), cutsites$Entry.name),]
cutsites_subset_unique_sequence_rank$gene_rank = cutsites_subset_unique_sequence_rank$sequence_rank

# extract the first instance of each unique gene from cutsites
# based on biological rank
cutsites = cutsites[order(cutsites$bio_rank),]
cutsites_subset_unique_bio_rank = cutsites[match(unique(cutsites$Gene_name), cutsites$Gene_name),]
cutsites_subset_unique_bio_rank$gene_rank = cutsites_subset_unique_bio_rank$bio_rank

# extract the first instance of each unique gene from cutsites
# based on fimo rank
cutsites = cutsites[order(cutsites$fimo_rank),]
cutsites_subset_unique_fimo_rank = cutsites[match(unique(cutsites$Gene_name), cutsites$Gene_name),]
cutsites_subset_unique_fimo_rank$gene_rank = cutsites_subset_unique_fimo_rank$fimo_rank

# extract GO2substrates rank columns of interest for visualisation
cutsites_subset_unique_GO2Substrates_extract = cutsites_subset_unique_GO2Substrates_rank %>% 
  select(Entry.name, Gene_name, gene_rank, literature_cleaved, tested_in_assay,
         cleaved_in_assay, Matched.sequence)
# append column to show rank type (to be used as a factor in spaghetti plot to follow)
cutsites_subset_unique_GO2Substrates_extract$rank_type = "GO2Substrates"
# append column for x axis order in spaghetti plot HARD CODED HERE
cutsites_subset_unique_GO2Substrates_extract$plot_order = 3
cutsites_subset_unique_GO2Substrates_extract$GO2Substrates_rank = cutsites_subset_unique_GO2Substrates_extract$gene_rank
cutsites_subset_unique_GO2Substrates_extract$fimo_rank = NA
cutsites_subset_unique_GO2Substrates_extract$bio_rank = NA
cutsites_subset_unique_GO2Substrates_extract$sequence_rank = NA

# extract sequence rank columns of interest for visualisation
cutsites_subset_unique_sequence_extract = cutsites_subset_unique_sequence_rank %>% 
  select(Entry.name, Gene_name, gene_rank, literature_cleaved, tested_in_assay,
         cleaved_in_assay, Matched.sequence)
# append column to show rank type (to be used as a factor in spaghetti plot to follow)
cutsites_subset_unique_sequence_extract$rank_type = "sequence"
# append column for x axis order in spaghetti plot HARD CODED HERE
cutsites_subset_unique_sequence_extract$plot_order = 2
cutsites_subset_unique_sequence_extract$GO2Substrates_rank = NA
cutsites_subset_unique_sequence_extract$fimo_rank = NA
cutsites_subset_unique_sequence_extract$bio_rank = NA
cutsites_subset_unique_sequence_extract$sequence_rank = cutsites_subset_unique_sequence_extract$gene_rank

# extract columns of interest for visualisation
cutsites_subset_unique_bio_extract = cutsites_subset_unique_bio_rank %>% 
  select(Entry.name, Gene_name, gene_rank, literature_cleaved, tested_in_assay,
         cleaved_in_assay, Matched.sequence)
# append column to show rank type (to be used as a factor in spaghetti plot to follow)
cutsites_subset_unique_bio_extract$rank_type = "bio"
cutsites_subset_unique_bio_extract$plot_order = 4
cutsites_subset_unique_bio_extract$GO2Substrates_rank = NA
cutsites_subset_unique_bio_extract$fimo_rank = NA
cutsites_subset_unique_bio_extract$bio_rank = cutsites_subset_unique_bio_extract$gene_rank
cutsites_subset_unique_bio_extract$sequence_rank = NA

# extract columns of interest for visualisation
cutsites_subset_unique_fimo_extract = cutsites_subset_unique_fimo_rank %>% 
  select(Entry.name, Gene_name, gene_rank, literature_cleaved, tested_in_assay,
         cleaved_in_assay, Matched.sequence)
# append column to show rank type (to be used as a factor in spaghetti plot to follow)
cutsites_subset_unique_fimo_extract$rank_type = "fimo"
cutsites_subset_unique_fimo_extract$plot_order = 1
cutsites_subset_unique_fimo_extract$GO2Substrates_rank = NA
cutsites_subset_unique_fimo_extract$fimo_rank = cutsites_subset_unique_fimo_extract$gene_rank
cutsites_subset_unique_fimo_extract$bio_rank = NA
cutsites_subset_unique_fimo_extract$sequence_rank = NA


# rbind final, sequence and bio dataframes to allow for spaghetti plot
spaghetti_cutsites = rbind.data.frame(cutsites_subset_unique_GO2Substrates_extract, cutsites_subset_unique_sequence_extract,
                                      cutsites_subset_unique_bio_extract, cutsites_subset_unique_fimo_extract)


# need to order extracted columns exactly the same, to allow calculation of differences
# between rank types

# loop: subset rows for each gene, populate new column with differences
# later, can remove rows with column of interest=0

for(i in 1:nrow(spaghetti_cutsites)) {
  if(!is.na(spaghetti_cutsites$Gene_name[i])) {
    k = spaghetti_cutsites %>% filter(Gene_name == Gene_name[i])
    spaghetti_cutsites$GO2Substrates_rank[i] = k$GO2Substrates_rank[which(k$rank_type=="GO2Substrates")] 
    spaghetti_cutsites$fimo_rank[i] = k$fimo_rank[which(k$rank_type=="fimo")] 
    spaghetti_cutsites$bio_rank[i] = k$bio_rank[which(k$rank_type=="bio")]
    spaghetti_cutsites$sequence_rank[i] = k$sequence_rank[which(k$rank_type=="sequence")] 
    print(i)
  }
}
# now extract table with single row for each, with all ranks
spaghetti_cutsites_unique = spaghetti_cutsites[match(unique(spaghetti_cutsites$Gene_name), 
                                                     spaghetti_cutsites$Gene_name),]

# extract columns of interest for visualisation
spaghetti_cutsites_unique_extract = spaghetti_cutsites_unique %>% 
  select(Gene_name, GO2Substrates_rank, bio_rank, sequence_rank)

cutsites_2 = cutsites

for(i in 1:nrow(cutsites_2)) {
  if(is.na(cutsites_2$cleaved_in_assay[i])) {
    cutsites_2$cleaved_in_assay[i] = "0"
  }
}

for(i in 1:nrow(cutsites_2)) {
  if(is.na(cutsites_2$tested_in_assay[i])) {
    cutsites_2$tested_in_assay[i] = "0"
  }
}

for(i in 1:nrow(cutsites_2)) {
  if(is.na(cutsites_2$literature_cleaved[i])) {
    cutsites_2$literature_cleaved[i] = "0"
  }
}

cutsites_2$sequence_cut_specific = NA
for(i in 1:nrow(cutsites_2)) {
  if(!is.na(cutsites$outcome_test_data[i])) {
    if(cutsites$outcome_test_data[i] == 1) {
      if(cutsites$sequence_cut[i] == as.character(cutsites$Matched.sequence[i])){
        cutsites_2$sequence_cut_specific[i] = "1"
      }
    }
  }
}

for(i in 1:nrow(cutsites_2)) {
  if(is.na(cutsites_2$sequence_cut_specific[i])) {
    cutsites_2$sequence_cut_specific[i] = "0"
  }
}


##### Make plots
# load required libraries
library(ggplot2)
library(gghighlight)
library(ggrepel)
library(lemon)

# Plot Fig. 3a
spaghetti_cutsites %>% filter(rank_type == "fimo" | 
                                rank_type == "bio") %>%
  ggplot(aes(x=reorder(rank_type, plot_order), gene_rank, 
             group = Gene_name, colour = Gene_name)) +
  scale_y_reverse(name = "Protein rank", breaks = c(0,500,1000,1500,2000,2500,3000,3500),
                  labels=c("0","500","1000","1500","2000","2500","3000","3500"),
                  limits=c(3500,0)) +
  coord_capped_cart(left='top', top=brackets_horizontal()) +
  theme(panel.grid.major.x = element_blank()) +
  scale_x_discrete(position = 'top') +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y=element_text(size = 8, face = "bold"))+
  theme(axis.text.x=element_text(size = 8, face = "bold")) +
  theme(axis.text.y=element_text(size = 8, face = "bold")) +
  theme(axis.line = element_line(color = "black",size = 0.5)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank()) +
  gghighlight(literature_cleaved == "1") +  
  geom_point(shape = "-", size = 12) +
  geom_line() +
  ggsave(width = 7, height = 12, units = "cm", filename ="Fig_3a.pdf")

# plot Fig. 3b
spaghetti_cutsites %>% filter(rank_type == "fimo" | 
                                rank_type == "sequence") %>%
  ggplot(aes(x=reorder(rank_type, plot_order), gene_rank, 
             group = Gene_name, colour = Gene_name)) +
  scale_y_reverse(name = "Protein rank", breaks = c(0,500,1000,1500,2000,2500,3000,3500),
                  labels=c("0","500","1000","1500","2000","2500","3000","3500"),
                  limits=c(3500,0)) +
  coord_capped_cart(left='top', top=brackets_horizontal()) +
  theme(panel.grid.major.x = element_blank()) +
  scale_x_discrete(position = 'top') +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y=element_text(size = 8, face = "bold"))+
  theme(axis.text.x=element_text(size = 8, face = "bold")) +
  theme(axis.text.y=element_text(size = 8, face = "bold")) +
  theme(axis.line = element_line(color = "black",size = 0.5)) +
  theme(legend.position = "none") +
  theme(panel.background = element_blank()) +
  gghighlight(literature_cleaved == "1") +  
  geom_point(shape = "-", size = 12) +
  geom_line() +
  ggsave(width =7, height = 12, units = "cm", filename = "Fig_3b.pdf")

# Plot Fig. 4a
spaghetti_cutsites %>% filter(rank_type == "GO2Substrates") %>%
  ggplot(aes(x=reorder(rank_type, plot_order), gene_rank, 
             group = Gene_name, colour = Gene_name)) +
  theme(legend.position = "none") +
  scale_y_reverse(name = "Protein rank", breaks = c(0,10,20,30),
                  labels=c("0","10","20","30"),
                  limits=c(30,0)) +
  coord_capped_cart(left='top', top=brackets_horizontal()) +
  theme(panel.grid.major.x = element_blank()) +
  scale_x_discrete(position = 'top') +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y=element_text(size = 8, face = "bold"))+
  theme(axis.text.x=element_text(size = 8, face = "bold")) +
  theme(axis.text.y=element_text(size = 8, face = "bold")) +
  theme(axis.line = element_line(color = "black",size = 0.5)) +
  theme(legend.text = element_text(size = 8, face = "bold")) +
  theme(panel.background = element_blank()) +
  gghighlight(literature_cleaved == "1") +
  geom_point(shape = "-", size = 28) +
  geom_line() +
  ggsave(width = 3, height = 12, units = "cm", filename ="Fig_4a.pdf")

### Plot Fig. 4b
sp=  ggplot(cutsites_2 %>%
              arrange(cleaved_in_assay),
            aes(sequence_rank_normalised, bio_rank_normalised, colour = tested_in_assay, alpha=tested_in_assay)) +
  scale_x_continuous(name = "Normalised Sequence Rank", position = 'top') +
  scale_y_reverse(name = "Normalised Function Rank") +
  theme(panel.background = element_blank()) +
  geom_point(size=2) +
  gghighlight(tested_in_assay=="1") +
  theme(legend.position = "bottom") +
  theme(legend.box.background = element_rect(colour = "black")) 
sp + scale_color_manual(values=c("red"))
ggsave(width = 12, height = 12, units = "cm", filename ="Fig_4b.pdf")


### Plot Fig. 4c (full plot)
# first apply unity-based normalisation to go-2-substrates scores
spaghetti_cutsites$GO2Substrates_rank_normalised = (spaghetti_cutsites$GO2Substrates_rank - min(spaghetti_cutsites$GO2Substrates_rank)) / (max(spaghetti_cutsites$GO2Substrates_rank) - min(spaghetti_cutsites$GO2Substrates_rank))

spaghetti_cutsites %>% filter(rank_type == "GO2Substrates") %>%
  ggplot(aes(x=reorder(rank_type, plot_order), GO2Substrates_rank_normalised, 
             group = Gene_name, colour = Gene_name)) +
  theme(legend.position = "none") +
  scale_y_reverse(name = "Protein rank", breaks = c(0,0.25,0.50, 0.75,1.00),
                  labels=c("0","0.25","0.5","0.75","1"),
                  limits=c(1,0)) +
  coord_capped_cart(left='top', top=brackets_horizontal()) +
  theme(panel.grid.major.x = element_blank()) +
  scale_x_discrete(position = 'top') +
  theme(panel.background = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y=element_text(size = 8, face = "bold"))+
  theme(axis.text.x=element_text(size = 8, face = "bold")) +
  theme(axis.text.y=element_text(size = 8, face = "bold")) +
  theme(axis.line = element_line(color = "black",size = 0.5)) +
  theme(legend.text = element_text(size = 8, face = "bold")) +
  theme(panel.background = element_blank()) +
  gghighlight(tested_in_assay == "1") +
  geom_point(fill = "red", size = 2, color = "red") +
  geom_point(fill = "black", size = 0.5, color = "red")
ggsave(width = 12, height = 12, units = "cm", filename ="Fig_4c_left.pdf")

# Plot Fig. 4c (zoom)
spaghetti_cutsites %>% filter(rank_type == "GO2Substrates") %>%
  ggplot(aes(x=reorder(rank_type, plot_order), GO2Substrates_rank_normalised, 
             group = Gene_name, colour = Gene_name)) +
  theme(legend.position = "none") +
  scale_y_reverse(name = "Protein rank", breaks = c(0,0.02,0.04, 0.06,0.08,0.1),
                  labels=c("0","0.02","0.04","0.06","0.08", "0.1"),
                  limits=c(0.1,0)) +
  coord_capped_cart(left='top', top=brackets_horizontal()) +
  theme(panel.grid.major.x = element_blank()) +
  scale_x_discrete(position = 'top') +
  theme(panel.background = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y=element_text(size = 8, face = "bold"))+
  theme(axis.text.x=element_text(size = 8, face = "bold")) +
  theme(axis.text.y=element_text(size = 8, face = "bold")) +
  theme(axis.line = element_line(color = "black",size = 0.5)) +
  theme(legend.text = element_text(size = 8, face = "bold")) +
  theme(panel.background = element_blank()) +
  gghighlight(tested_in_assay == "1") +
  geom_point(fill = "red", size = 2, color = "red") +
  geom_point(fill = "black", size = 0.5, color = "red")
ggsave(width = 12, height = 12, units = "cm", filename ="Fig_4c_right.pdf")

### Generate ROC and PR curves and plot for Fig 8c-e and h-k
library(precrec)

#calculate ROC and precision-recall for function (bio) module
subset_bio = cutsites_subset_unique_GO2Substrates_rank %>% filter(!is.na(cutsites_subset_unique_GO2Substrates_rank$outcome_test_data))
subset_bio = subset(cutsites_subset_unique_GO2Substrates_rank, select = c(outcome_test_data, bio_rank_normalised))

#note - have to reverse normalised rank here (1= best) precrec doesnt allow direction setting
subset_bio$bio_rank_normalised = 1-subset_bio$bio_rank_normalised

tested_indices_bio = which(!is.na(subset_bio$outcome_test_data))
msmdat_bio <- mmdata(subset_bio$bio_rank_normalised[tested_indices_bio], 
                     subset_bio$outcome_test_data[tested_indices_bio])
sscurves <- evalmod(msmdat_bio)
# plot Fig 8c and 8h
pdf("Fig_8c_8h.pdf") 
# 2. Create a plot
plot(sscurves)
# Close the pdf file
dev.off()


#calculate ROC and precision-recall for sequence module 
subset_sequence = cutsites_subset_unique_GO2Substrates_rank %>% filter(!is.na(cutsites_subset_unique_GO2Substrates_rank$outcome_test_data))
subset_sequence = subset(cutsites_subset_unique_GO2Substrates_rank, select = c(outcome_test_data, sequence_rank_normalised))

#note - have to reverse normalised rank here (1= best) precrec doesnt allow direction setting
subset_sequence$sequence_rank_normalised = 1-subset_sequence$sequence_rank_normalised

tested_indices_sequence = which(!is.na(subset_sequence$outcome_test_data))
msmdat_seq <- mmdata(subset_sequence$sequence_rank_normalised[tested_indices_sequence], 
                     subset_sequence$outcome_test_data[tested_indices_sequence])
sscurves <- evalmod(msmdat_seq)
# plot Fig 8d and 8i
pdf("Fig_8d_8i.pdf") 
# 2. Create a plot
plot(sscurves)
# Close the pdf file
dev.off()


#calculate ROC and precision-recall for go2substrates
subset_go2substrates = cutsites_subset_unique_GO2Substrates_rank %>% filter(!is.na(cutsites_subset_unique_GO2Substrates_rank$outcome_test_data))
subset_go2substrates = subset(cutsites_subset_unique_GO2Substrates_rank, select = c(outcome_test_data, final_rank_normalised))

#note - have to reverse normalised rank here (1= best) precrec doesnt allow direction setting
subset_go2substrates$final_rank_normalised = 1-subset_go2substrates$final_rank_normalised
tested_indices_go2substrates = which(!is.na(subset_go2substrates$outcome_test_data))
msmdat_go2substrates <- mmdata(subset_go2substrates$final_rank_normalised[tested_indices_go2substrates], 
                               subset_go2substrates$outcome_test_data[tested_indices_go2substrates])
sscurves <- evalmod(msmdat_go2substrates)
# plot Fig 8e and 8h
pdf("Fig_8e_8j.pdf") 
# 2. Create a plot
plot(sscurves)
# Close the pdf file
dev.off()



#### Calculate precision vs rank info for Fig 8i, and plot Fig 8k #####
test_subset = cutsites_subset_unique_GO2Substrates_rank %>% filter(!is.na(cutsites_subset_unique_GO2Substrates_rank$outcome_test_data))
test_subset2 = subset(test_subset, select = c(cleaved_in_assay, outcome_test_data, gene_rank))
test_subset2 = test_subset2[order(test_subset2$gene_rank, decreasing = FALSE),]

test_subset2$TP = NA
test_subset2 = rbind(data.frame(cleaved_in_assay = NA, outcome_test_data = NA, gene_rank = 0, TP=0), test_subset2)

# identify first row where outcome_test_data = 1
temp_constant = min(which(test_subset2$outcome_test_data == 1))

# calculate number of true positives at each rank tested
for(i in temp_constant:nrow(test_subset2)) {
  if(test_subset2$outcome_test_data[i] == 1){
    test_subset2$TP[i] = max(test_subset2$TP, na.rm=TRUE) + 1
  }
}

for(i in temp_constant:nrow(test_subset2)) {
  if(is.na(test_subset2$TP[i])){
    test_subset2$TP[i] = max(test_subset2$TP[1:i], na.rm = TRUE)
  }
}

# calculate number of false positives at each rank tested
test_subset2$FP = NA
test_subset2$FP[1] = 0
# identify first row where outcome_test_data = 1
temp_constant = min(which(test_subset2$outcome_test_data == 1))
# calculate number of false positives at each rank tested
for(i in 2:nrow(test_subset2)) {
  if(test_subset2$outcome_test_data[i] == 0){
    test_subset2$FP[i] = max(test_subset2$FP, na.rm=TRUE) + 1
  }
}

for(i in temp_constant:nrow(test_subset2)) {
  if(is.na(test_subset2$FP[i])){
    test_subset2$FP[i] = max(test_subset2$FP[1:i], na.rm = TRUE)
  }
}

# calculate number of true negatives at each rank tested
test_subset2$TN = NA
test_subset2$TN[1] = 0
# calculate number of true negatives at each rank tested
for(i in 2:nrow(test_subset2)) {
  #if(test_subset2$outcome_test_data[i] == 1){
  test_subset2$TN[i] = length(which(test_subset2$outcome_test_data == 0)) - test_subset2$FP[i]
  #length(which(test_subset2$outcome_test_data == 0)) - test_subset2$TP[i]
}
#}

# calculate number of false negatives at each rank tested
test_subset2$FN = NA
test_subset2$FN[1] = max(test_subset2$TP)

for(i in temp_constant:nrow(test_subset2)) {
  if(test_subset2$outcome_test_data[i] == 1){
    test_subset2$FN[i] = min(test_subset2$FN, na.rm=TRUE) - 1
  }
}

for(i in temp_constant:nrow(test_subset2)) {
  if(is.na(test_subset2$FN[i])){
    test_subset2$FN[i] = min(test_subset2$FN[1:i], na.rm = TRUE)
  }
}

# calculate precision at each rank
test_subset2$precision = test_subset2$TP/(test_subset2$TP+test_subset2$FP)
test_subset2$sensitivity = test_subset2$TP/(test_subset2$TP + test_subset2$FN)
test_subset2$specificity = test_subset2$TN/(test_subset2$TN + test_subset2$FP)

# setup $cleaved_in_assay as a logical to allow compatibility with scale_color_manual
test_subset2$cleaved_in_assay = as.logical(test_subset2$cleaved_in_assay)
test_subset2[is.na(test_subset2)] = FALSE
# remove top row, since it now has a value and will appear as 0,0 on plot
test_subset2 = test_subset2[-1,]

### Plot Fig. 8k (full plot)
test_subset2 %>% ggplot(aes(gene_rank, precision,colour = cleaved_in_assay)) +
  scale_y_continuous(limits = c(0,1.02), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0,4000), expand = c(0, 0),
                     name = "GO-2-Substrates Rank", breaks = c(0,500,1000,1500,2000,2500,3000,3500,4000),
                     labels=c("0","500","1000","1500","2000","2500","3000","3500","4000")) +
  geom_point(size=2)+
  scale_color_manual(values=c("black", "red")) +
  geom_hline(yintercept = 0.2, linetype = 2) +
  theme(panel.background = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))
ggsave(width = 15, height = 10, units = "cm", filename ="Fig_8k_full.pdf")

### Plot Fig. 8k (subset plot in paper)
test_subset2 %>% ggplot(aes(gene_rank, precision,colour = cleaved_in_assay)) +
  scale_y_continuous(limits = c(0,1.02), expand = c(0, 0), 
                     name = "Precision") +
  scale_x_continuous(limits = c(0,200), expand = c(0, 0),
                     name = "GO-2-Substrates Rank", breaks = c(0,20,40,60,80,100,120,140,160,180,200),
                     labels=c("0","20","40","60","80","100","120","140","160","180","200")) +
  geom_point(size=2)+
  scale_color_manual(values=c("black", "red")) +
  geom_hline(yintercept = 0.2, linetype = 2) +
  theme(panel.background = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))
ggsave(width = 15, height = 10, units = "cm", filename ="Fig_8k.pdf")


###### Plot Supplementary Fig 1 #########
### note - geom_abline hard coded here with optimal linear model derived in script A
# as follows: fit_for_plot = -5.460511e-05(x) + 0.7760202
test_plot=  ggplot(ordered_GO_BP_string, aes(background.gene.count, strength)) +
  scale_x_continuous(name = "Number of proteins annotated with GOTERM", limits = c(0, 10000), 
                     expand = c(0,0)) +
  scale_y_continuous(name = "Strength of enrichment", limits=c(0,3),
                     expand = c(0,0), 
                     breaks = c(0,0.5,1,1.5,2,2.5,3),
                     labels=c("0","0.5","1","1.5","2","2.5","3"),) +
  theme(panel.background = element_blank()) +
  theme(axis.line.x = element_line(colour="black", size = 0.5)) +
  theme(axis.line.y = element_line(colour="black", size = 0.5)) +
  geom_point(size=1, colour = "steelblue") +
  geom_abline(slope = -5.460511e-05, intercept = 0.7760202, colour="red") +
  theme(legend.position = "bottom") +
  theme(legend.box.background = element_rect(colour = "black")) +
  geom_hline(yintercept = 1.43, linetype = 2) +
  geom_hline(yintercept = 0.7760202, linetype = 2) 
test_plot
ggsave(width = 15, height = 10, units = "cm", filename ="Supplementary_Fig.1.pdf")


###### Prepare Data for Supplementary Fig 2 #######
# calculate GOTERMs belonging to each criteria
new_GO_BP_string_criteria_1 = new_GO_BP_string[new_GO_BP_string$X.term.ID %in% GO_BP_shell1,]
new_GO_BP_string_criteria_2 = new_GO_BP_string[new_GO_BP_string$X.term.ID %in% GO_BP_shell2,]
temp_value = new_GO_BP_string_criteria_2$X.term.ID %!in% GO_BP_shell1
new_GO_BP_string_criteria_2 = new_GO_BP_string_criteria_2[temp_value,]
new_GO_BP_string_criteria_3 = new_GO_BP_string[new_GO_BP_string$X.term.ID %in% GO_BP_shell3,]
temp_value = new_GO_BP_string_criteria_3$X.term.ID %!in% c(GO_BP_shell1, new_GO_BP_string_criteria_2$X.term.ID)
new_GO_BP_string_criteria_3 = new_GO_BP_string_criteria_3[temp_value,]

# calculate for criteria 1
TNFAIP3_criteria_1 = sum(str_count(new_GO_BP_string_criteria_1$matching.proteins.in.your.network..labels., "TNFAIP3"))
BCL10_criteria_1 = sum(str_count(new_GO_BP_string_criteria_1$matching.proteins.in.your.network..labels., "BCL10"))
N4BP1_criteria_1 = sum(str_count(new_GO_BP_string_criteria_1$matching.proteins.in.your.network..labels., "N4BP1"))
RC3H1_criteria_1 = sum(str_count(new_GO_BP_string_criteria_1$matching.proteins.in.your.network..labels., "RC3H1"))
RC3H2_criteria_1 = sum(str_count(new_GO_BP_string_criteria_1$matching.proteins.in.your.network..labels., "RC3H2"))
MALT1_criteria_1 = sum(str_count(new_GO_BP_string_criteria_1$matching.proteins.in.your.network..labels., "MALT1"))
LIMA1_criteria_1 = sum(str_count(new_GO_BP_string_criteria_1$matching.proteins.in.your.network..labels., "LIMA1"))
MAP3K14_criteria_1 = sum(str_count(new_GO_BP_string_criteria_1$matching.proteins.in.your.network..labels., "MAP3K14"))
RELB_criteria_1 = sum(str_count(new_GO_BP_string_criteria_1$matching.proteins.in.your.network..labels., "RELB"))
RBCK1_criteria_1 = sum(str_count(new_GO_BP_string_criteria_1$matching.proteins.in.your.network..labels., "RBCK1"))
CYLD_criteria_1 = sum(str_count(new_GO_BP_string_criteria_1$matching.proteins.in.your.network..labels., "CYLD"))
ZC3H12A_criteria_1 = sum(str_count(new_GO_BP_string_criteria_1$matching.proteins.in.your.network..labels., "ZC3H12A"))

# calculate for criteria 2
TNFAIP3_criteria_2 = sum(str_count(new_GO_BP_string_criteria_2$matching.proteins.in.your.network..labels., "TNFAIP3"))
BCL10_criteria_2 = sum(str_count(new_GO_BP_string_criteria_2$matching.proteins.in.your.network..labels., "BCL10"))
N4BP1_criteria_2 = sum(str_count(new_GO_BP_string_criteria_2$matching.proteins.in.your.network..labels., "N4BP1"))
RC3H1_criteria_2 = sum(str_count(new_GO_BP_string_criteria_2$matching.proteins.in.your.network..labels., "RC3H1"))
RC3H2_criteria_2 = sum(str_count(new_GO_BP_string_criteria_2$matching.proteins.in.your.network..labels., "RC3H2"))
MALT1_criteria_2 = sum(str_count(new_GO_BP_string_criteria_2$matching.proteins.in.your.network..labels., "MALT1"))
LIMA1_criteria_2 = sum(str_count(new_GO_BP_string_criteria_2$matching.proteins.in.your.network..labels., "LIMA1"))
MAP3K14_criteria_2 = sum(str_count(new_GO_BP_string_criteria_2$matching.proteins.in.your.network..labels., "MAP3K14"))
RELB_criteria_2 = sum(str_count(new_GO_BP_string_criteria_2$matching.proteins.in.your.network..labels., "RELB"))
RBCK1_criteria_2 = sum(str_count(new_GO_BP_string_criteria_2$matching.proteins.in.your.network..labels., "RBCK1"))
CYLD_criteria_2 = sum(str_count(new_GO_BP_string_criteria_2$matching.proteins.in.your.network..labels., "CYLD"))
ZC3H12A_criteria_2 = sum(str_count(new_GO_BP_string_criteria_2$matching.proteins.in.your.network..labels., "ZC3H12A"))


# calculate for criteria 3
TNFAIP3_criteria_3 = sum(str_count(new_GO_BP_string_criteria_3$matching.proteins.in.your.network..labels., "TNFAIP3"))
BCL10_criteria_3 = sum(str_count(new_GO_BP_string_criteria_3$matching.proteins.in.your.network..labels., "BCL10"))
N4BP1_criteria_3 = sum(str_count(new_GO_BP_string_criteria_3$matching.proteins.in.your.network..labels., "N4BP1"))
RC3H1_criteria_3 = sum(str_count(new_GO_BP_string_criteria_3$matching.proteins.in.your.network..labels., "RC3H1"))
RC3H2_criteria_3 = sum(str_count(new_GO_BP_string_criteria_3$matching.proteins.in.your.network..labels., "RC3H2"))
MALT1_criteria_3 = sum(str_count(new_GO_BP_string_criteria_3$matching.proteins.in.your.network..labels., "MALT1"))
LIMA1_criteria_3 = sum(str_count(new_GO_BP_string_criteria_3$matching.proteins.in.your.network..labels., "LIMA1"))
MAP3K14_criteria_3 = sum(str_count(new_GO_BP_string_criteria_3$matching.proteins.in.your.network..labels., "MAP3K14"))
RELB_criteria_3 = sum(str_count(new_GO_BP_string_criteria_3$matching.proteins.in.your.network..labels., "RELB"))
RBCK1_criteria_3 = sum(str_count(new_GO_BP_string_criteria_3$matching.proteins.in.your.network..labels., "RBCK1"))
CYLD_criteria_3 = sum(str_count(new_GO_BP_string_criteria_3$matching.proteins.in.your.network..labels., "CYLD"))
ZC3H12A_criteria_3 = sum(str_count(new_GO_BP_string_criteria_3$matching.proteins.in.your.network..labels., "ZC3H12A"))


# assemble data frame for plotting with ggplot2 
# note... very hard coded! need to be quick

# first, initialise data frame data into proteins
BP_data_frame = data.frame(NA, NA, NA)
BP_data_frame [ nrow(BP_data_frame) + 35 , ] <- NA
colnames(BP_data_frame) = c("Protein", "Number of annotations", "Criteria")

#set criteria number for each
for(i in seq(1, nrow(BP_data_frame)-2, 3)) {
  BP_data_frame$Criteria[i] = 1
  BP_data_frame$Criteria[i+1] = 2
  BP_data_frame$Criteria[i+2] = 3
}

# set protein name for each
Protein_names = c("TNFAIP3", "BCL10", "N4BP1", "RC3H1", "RC3H2", "MALT1", "LIMA1", "MAP3K14", "RELB", "RBCK1", "CYLD", "ZC3H12A")
for(i in seq(4, nrow(BP_data_frame), 3)) {
  constant = (i + 2) / 3  
  BP_data_frame$Protein[i] = Protein_names[constant]
  BP_data_frame$Protein[i+1] = Protein_names[constant]
  BP_data_frame$Protein[i+2] = Protein_names[constant]
}

# now addnumber of annotations calculated previously
annotations = c(TNFAIP3_criteria_1, TNFAIP3_criteria_2, TNFAIP3_criteria_3, 
                BCL10_criteria_1, BCL10_criteria_2, BCL10_criteria_3, 
                N4BP1_criteria_1, N4BP1_criteria_2, N4BP1_criteria_3, 
                RC3H1_criteria_1, RC3H1_criteria_2, RC3H1_criteria_3,
                RC3H2_criteria_1, RC3H2_criteria_2, RC3H2_criteria_3, 
                MALT1_criteria_1, MALT1_criteria_2, MALT1_criteria_3, 
                LIMA1_criteria_1, LIMA1_criteria_2, LIMA1_criteria_3, 
                MAP3K14_criteria_1, MAP3K14_criteria_2, MAP3K14_criteria_3,
                RELB_criteria_1, RELB_criteria_2, RELB_criteria_3, 
                RBCK1_criteria_1, RBCK1_criteria_2, RBCK1_criteria_3, 
                CYLD_criteria_1, CYLD_criteria_2, CYLD_criteria_3, 
                ZC3H12A_criteria_1, ZC3H12A_criteria_2, ZC3H12A_criteria_3)

# add info to dataframe
BP_data_frame[2] = annotations 
colnames(BP_data_frame)[2] = "Number_of_annotations"
BP_data_frame[2] = as.numeric(unlist(BP_data_frame[2]))
BP_data_frame[3] = as.character(unlist(BP_data_frame[3]))
# now plot with ggplot2
library(ggplot2)
ggplot(BP_data_frame, aes(x=Protein, y=Number_of_annotations, fill=Criteria)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_y_continuous(limits = c(0,60), breaks = c(0, 10, 20, 30, 40, 50, 60), labels = c(0, 10, 20, 30, 40, 50, 60)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x=element_text(size = 12, face = "bold")) +
  theme(axis.title.y=element_text(size = 12, face = "bold")) +
  theme(axis.text.x=element_text(size = 12, face = "bold")) +
  theme(axis.text.y=element_text(size = 12, face = "bold")) +
  theme(axis.line.y = element_line(color = "black",size = 0.5)) +
  theme(legend.text = element_text(size = 12, face = "bold")) +
  theme(legend.title = element_text(size = 12, face = "bold")) + 
  theme(axis.ticks.x = element_blank())
ggsave(width = 30, height = 7, units = "cm", filename ="Supplementary_Fig2_GO_BP.pdf")


# repeat for other GO term types

# calculate GOTERMs belonging to each criteria
new_GO_MF_string_criteria_1 = new_GO_MF_string[new_GO_MF_string$X.term.ID %in% GO_MF_shell1,]
new_GO_MF_string_criteria_2 = new_GO_MF_string[new_GO_MF_string$X.term.ID %in% GO_MF_shell2,]
temp_value = new_GO_MF_string_criteria_2$X.term.ID %!in% GO_MF_shell1
new_GO_MF_string_criteria_2 = new_GO_MF_string_criteria_2[temp_value,]
new_GO_MF_string_criteria_3 = new_GO_MF_string[new_GO_MF_string$X.term.ID %in% GO_MF_shell3,]
temp_value = new_GO_MF_string_criteria_3$X.term.ID %!in% c(GO_MF_shell1, new_GO_MF_string_criteria_2$X.term.ID)
new_GO_MF_string_criteria_3 = new_GO_MF_string_criteria_3[temp_value,]

# calculate for criteria 1
TNFAIP3_criteria_1 = sum(str_count(new_GO_MF_string_criteria_1$matching.proteins.in.your.network..labels., "TNFAIP3"))
BCL10_criteria_1 = sum(str_count(new_GO_MF_string_criteria_1$matching.proteins.in.your.network..labels., "BCL10"))
N4BP1_criteria_1 = sum(str_count(new_GO_MF_string_criteria_1$matching.proteins.in.your.network..labels., "N4BP1"))
RC3H1_criteria_1 = sum(str_count(new_GO_MF_string_criteria_1$matching.proteins.in.your.network..labels., "RC3H1"))
RC3H2_criteria_1 = sum(str_count(new_GO_MF_string_criteria_1$matching.proteins.in.your.network..labels., "RC3H2"))
MALT1_criteria_1 = sum(str_count(new_GO_MF_string_criteria_1$matching.proteins.in.your.network..labels., "MALT1"))
LIMA1_criteria_1 = sum(str_count(new_GO_MF_string_criteria_1$matching.proteins.in.your.network..labels., "LIMA1"))
MAP3K14_criteria_1 = sum(str_count(new_GO_MF_string_criteria_1$matching.proteins.in.your.network..labels., "MAP3K14"))
RELB_criteria_1 = sum(str_count(new_GO_MF_string_criteria_1$matching.proteins.in.your.network..labels., "RELB"))
RBCK1_criteria_1 = sum(str_count(new_GO_MF_string_criteria_1$matching.proteins.in.your.network..labels., "RBCK1"))
CYLD_criteria_1 = sum(str_count(new_GO_MF_string_criteria_1$matching.proteins.in.your.network..labels., "CYLD"))
ZC3H12A_criteria_1 = sum(str_count(new_GO_MF_string_criteria_1$matching.proteins.in.your.network..labels., "ZC3H12A"))

# calculate for criteria 2
TNFAIP3_criteria_2 = sum(str_count(new_GO_MF_string_criteria_2$matching.proteins.in.your.network..labels., "TNFAIP3"))
BCL10_criteria_2 = sum(str_count(new_GO_MF_string_criteria_2$matching.proteins.in.your.network..labels., "BCL10"))
N4BP1_criteria_2 = sum(str_count(new_GO_MF_string_criteria_2$matching.proteins.in.your.network..labels., "N4BP1"))
RC3H1_criteria_2 = sum(str_count(new_GO_MF_string_criteria_2$matching.proteins.in.your.network..labels., "RC3H1"))
RC3H2_criteria_2 = sum(str_count(new_GO_MF_string_criteria_2$matching.proteins.in.your.network..labels., "RC3H2"))
MALT1_criteria_2 = sum(str_count(new_GO_MF_string_criteria_2$matching.proteins.in.your.network..labels., "MALT1"))
LIMA1_criteria_2 = sum(str_count(new_GO_MF_string_criteria_2$matching.proteins.in.your.network..labels., "LIMA1"))
MAP3K14_criteria_2 = sum(str_count(new_GO_MF_string_criteria_2$matching.proteins.in.your.network..labels., "MAP3K14"))
RELB_criteria_2 = sum(str_count(new_GO_MF_string_criteria_2$matching.proteins.in.your.network..labels., "RELB"))
RBCK1_criteria_2 = sum(str_count(new_GO_MF_string_criteria_2$matching.proteins.in.your.network..labels., "RBCK1"))
CYLD_criteria_2 = sum(str_count(new_GO_MF_string_criteria_2$matching.proteins.in.your.network..labels., "CYLD"))
ZC3H12A_criteria_2 = sum(str_count(new_GO_MF_string_criteria_2$matching.proteins.in.your.network..labels., "ZC3H12A"))


# calculate for criteria 3
TNFAIP3_criteria_3 = sum(str_count(new_GO_MF_string_criteria_3$matching.proteins.in.your.network..labels., "TNFAIP3"))
BCL10_criteria_3 = sum(str_count(new_GO_MF_string_criteria_3$matching.proteins.in.your.network..labels., "BCL10"))
N4BP1_criteria_3 = sum(str_count(new_GO_MF_string_criteria_3$matching.proteins.in.your.network..labels., "N4BP1"))
RC3H1_criteria_3 = sum(str_count(new_GO_MF_string_criteria_3$matching.proteins.in.your.network..labels., "RC3H1"))
RC3H2_criteria_3 = sum(str_count(new_GO_MF_string_criteria_3$matching.proteins.in.your.network..labels., "RC3H2"))
MALT1_criteria_3 = sum(str_count(new_GO_MF_string_criteria_3$matching.proteins.in.your.network..labels., "MALT1"))
LIMA1_criteria_3 = sum(str_count(new_GO_MF_string_criteria_3$matching.proteins.in.your.network..labels., "LIMA1"))
MAP3K14_criteria_3 = sum(str_count(new_GO_MF_string_criteria_3$matching.proteins.in.your.network..labels., "MAP3K14"))
RELB_criteria_3 = sum(str_count(new_GO_MF_string_criteria_3$matching.proteins.in.your.network..labels., "RELB"))
RBCK1_criteria_3 = sum(str_count(new_GO_MF_string_criteria_3$matching.proteins.in.your.network..labels., "RBCK1"))
CYLD_criteria_3 = sum(str_count(new_GO_MF_string_criteria_3$matching.proteins.in.your.network..labels., "CYLD"))
ZC3H12A_criteria_3 = sum(str_count(new_GO_MF_string_criteria_3$matching.proteins.in.your.network..labels., "ZC3H12A"))


# assemble data frame for plotting with ggplot2 
# note... very hard coded! need to be quick

# first, initialise data frame data into proteins
MF_data_frame = data.frame(NA, NA, NA)
MF_data_frame [ nrow(MF_data_frame) + 35 , ] <- NA
colnames(MF_data_frame) = c("Protein", "Number of annotations", "Criteria")
#set criteria number for each
for(i in seq(1, nrow(MF_data_frame)-2, 3)) {
  MF_data_frame$Criteria[i] = 1
  MF_data_frame$Criteria[i+1] = 2
  MF_data_frame$Criteria[i+2] = 3
}
# set protein name for each
Protein_names = c("TNFAIP3", "BCL10", "N4BP1", "RC3H1", "RC3H2", "MALT1", "LIMA1", "MAP3K14", "RELB", "RBCK1", "CYLD", "ZC3H12A")
for(i in seq(1, nrow(MF_data_frame), 3)) {
  constant = (i + 2) / 3  
  MF_data_frame$Protein[i] = Protein_names[constant]
  MF_data_frame$Protein[i+1] = Protein_names[constant]
  MF_data_frame$Protein[i+2] = Protein_names[constant]
}
# now addnumber of annotations calculated previously
annotations = c(TNFAIP3_criteria_1, TNFAIP3_criteria_2, TNFAIP3_criteria_3, 
                BCL10_criteria_1, BCL10_criteria_2, BCL10_criteria_3, 
                N4BP1_criteria_1, N4BP1_criteria_2, N4BP1_criteria_3, 
                RC3H1_criteria_1, RC3H1_criteria_2, RC3H1_criteria_3,
                RC3H2_criteria_1, RC3H2_criteria_2, RC3H2_criteria_3, 
                MALT1_criteria_1, MALT1_criteria_2, MALT1_criteria_3, 
                LIMA1_criteria_1, LIMA1_criteria_2, LIMA1_criteria_3, 
                MAP3K14_criteria_1, MAP3K14_criteria_2, MAP3K14_criteria_3,
                RELB_criteria_1, RELB_criteria_2, RELB_criteria_3, 
                RBCK1_criteria_1, RBCK1_criteria_2, RBCK1_criteria_3, 
                CYLD_criteria_1, CYLD_criteria_2, CYLD_criteria_3, 
                ZC3H12A_criteria_1, ZC3H12A_criteria_2, ZC3H12A_criteria_3)
# add info to dataframe
MF_data_frame[2] = annotations 
colnames(MF_data_frame)[2] = "Number_of_annotations"
MF_data_frame[2] = as.numeric(unlist(MF_data_frame[2]))
MF_data_frame[3] = as.character(unlist(MF_data_frame[3]))
# now plot with ggplot2
library(ggplot2)
ggplot(MF_data_frame, aes(x=Protein, y=Number_of_annotations, fill=Criteria)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_y_continuous(limits = c(0,6), breaks = c(0, 2, 4, 6), labels = c(0, 2, 4, 6)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x=element_text(size = 12, face = "bold")) +
  theme(axis.title.y=element_text(size = 12, face = "bold")) +
  theme(axis.text.x=element_text(size = 12, face = "bold")) +
  theme(axis.text.y=element_text(size = 12, face = "bold")) +
  theme(axis.line.y = element_line(color = "black",size = 0.5)) +
  theme(legend.text = element_text(size = 12, face = "bold")) +
  theme(legend.title = element_text(size = 12, face = "bold")) + 
  theme(axis.ticks.x = element_blank())
ggsave(width = 30, height = 7, units = "cm", filename ="Supplementary_Fig2_GO_MF.pdf")

#### repeat for GO_CC
# calculate GOTERMs belonging to each criteria
new_GO_CC_string_criteria_1 = new_GO_CC_string[new_GO_CC_string$X.term.ID %in% GO_CC_shell1,]
new_GO_CC_string_criteria_2 = new_GO_CC_string[new_GO_CC_string$X.term.ID %in% GO_CC_shell2,]
temp_value = new_GO_CC_string_criteria_2$X.term.ID %!in% GO_CC_shell1
new_GO_CC_string_criteria_2 = new_GO_CC_string_criteria_2[temp_value,]
new_GO_CC_string_criteria_3 = new_GO_CC_string[new_GO_CC_string$X.term.ID %in% GO_CC_shell3,]
temp_value = new_GO_CC_string_criteria_3$X.term.ID %!in% c(GO_CC_shell1, new_GO_CC_string_criteria_2$X.term.ID)
new_GO_CC_string_criteria_3 = new_GO_CC_string_criteria_3[temp_value,]

# calculate for criteria 1
TNFAIP3_criteria_1 = sum(str_count(new_GO_CC_string_criteria_1$matching.proteins.in.your.network..labels., "TNFAIP3"))
BCL10_criteria_1 = sum(str_count(new_GO_CC_string_criteria_1$matching.proteins.in.your.network..labels., "BCL10"))
N4BP1_criteria_1 = sum(str_count(new_GO_CC_string_criteria_1$matching.proteins.in.your.network..labels., "N4BP1"))
RC3H1_criteria_1 = sum(str_count(new_GO_CC_string_criteria_1$matching.proteins.in.your.network..labels., "RC3H1"))
RC3H2_criteria_1 = sum(str_count(new_GO_CC_string_criteria_1$matching.proteins.in.your.network..labels., "RC3H2"))
MALT1_criteria_1 = sum(str_count(new_GO_CC_string_criteria_1$matching.proteins.in.your.network..labels., "MALT1"))
LIMA1_criteria_1 = sum(str_count(new_GO_CC_string_criteria_1$matching.proteins.in.your.network..labels., "LIMA1"))
MAP3K14_criteria_1 = sum(str_count(new_GO_CC_string_criteria_1$matching.proteins.in.your.network..labels., "MAP3K14"))
RELB_criteria_1 = sum(str_count(new_GO_CC_string_criteria_1$matching.proteins.in.your.network..labels., "RELB"))
RBCK1_criteria_1 = sum(str_count(new_GO_CC_string_criteria_1$matching.proteins.in.your.network..labels., "RBCK1"))
CYLD_criteria_1 = sum(str_count(new_GO_CC_string_criteria_1$matching.proteins.in.your.network..labels., "CYLD"))
ZC3H12A_criteria_1 = sum(str_count(new_GO_CC_string_criteria_1$matching.proteins.in.your.network..labels., "ZC3H12A"))

# calculate for criteria 2
TNFAIP3_criteria_2 = sum(str_count(new_GO_CC_string_criteria_2$matching.proteins.in.your.network..labels., "TNFAIP3"))
BCL10_criteria_2 = sum(str_count(new_GO_CC_string_criteria_2$matching.proteins.in.your.network..labels., "BCL10"))
N4BP1_criteria_2 = sum(str_count(new_GO_CC_string_criteria_2$matching.proteins.in.your.network..labels., "N4BP1"))
RC3H1_criteria_2 = sum(str_count(new_GO_CC_string_criteria_2$matching.proteins.in.your.network..labels., "RC3H1"))
RC3H2_criteria_2 = sum(str_count(new_GO_CC_string_criteria_2$matching.proteins.in.your.network..labels., "RC3H2"))
MALT1_criteria_2 = sum(str_count(new_GO_CC_string_criteria_2$matching.proteins.in.your.network..labels., "MALT1"))
LIMA1_criteria_2 = sum(str_count(new_GO_CC_string_criteria_2$matching.proteins.in.your.network..labels., "LIMA1"))
MAP3K14_criteria_2 = sum(str_count(new_GO_CC_string_criteria_2$matching.proteins.in.your.network..labels., "MAP3K14"))
RELB_criteria_2 = sum(str_count(new_GO_CC_string_criteria_2$matching.proteins.in.your.network..labels., "RELB"))
RBCK1_criteria_2 = sum(str_count(new_GO_CC_string_criteria_2$matching.proteins.in.your.network..labels., "RBCK1"))
CYLD_criteria_2 = sum(str_count(new_GO_CC_string_criteria_2$matching.proteins.in.your.network..labels., "CYLD"))
ZC3H12A_criteria_2 = sum(str_count(new_GO_CC_string_criteria_2$matching.proteins.in.your.network..labels., "ZC3H12A"))


# calculate for criteria 3
TNFAIP3_criteria_3 = sum(str_count(new_GO_CC_string_criteria_3$matching.proteins.in.your.network..labels., "TNFAIP3"))
BCL10_criteria_3 = sum(str_count(new_GO_CC_string_criteria_3$matching.proteins.in.your.network..labels., "BCL10"))
N4BP1_criteria_3 = sum(str_count(new_GO_CC_string_criteria_3$matching.proteins.in.your.network..labels., "N4BP1"))
RC3H1_criteria_3 = sum(str_count(new_GO_CC_string_criteria_3$matching.proteins.in.your.network..labels., "RC3H1"))
RC3H2_criteria_3 = sum(str_count(new_GO_CC_string_criteria_3$matching.proteins.in.your.network..labels., "RC3H2"))
MALT1_criteria_3 = sum(str_count(new_GO_CC_string_criteria_3$matching.proteins.in.your.network..labels., "MALT1"))
LIMA1_criteria_3 = sum(str_count(new_GO_CC_string_criteria_3$matching.proteins.in.your.network..labels., "LIMA1"))
MAP3K14_criteria_3 = sum(str_count(new_GO_CC_string_criteria_3$matching.proteins.in.your.network..labels., "MAP3K14"))
RELB_criteria_3 = sum(str_count(new_GO_CC_string_criteria_3$matching.proteins.in.your.network..labels., "RELB"))
RBCK1_criteria_3 = sum(str_count(new_GO_CC_string_criteria_3$matching.proteins.in.your.network..labels., "RBCK1"))
CYLD_criteria_3 = sum(str_count(new_GO_CC_string_criteria_3$matching.proteins.in.your.network..labels., "CYLD"))
ZC3H12A_criteria_3 = sum(str_count(new_GO_CC_string_criteria_3$matching.proteins.in.your.network..labels., "ZC3H12A"))


# assemble data frame for plotting with ggplot2 
# note... very hard coded! need to be quick

# first, initialise data frame data into proteins
CC_data_frame = data.frame(NA, NA, NA)
CC_data_frame [ nrow(CC_data_frame) + 35 , ] <- NA
colnames(CC_data_frame) = c("Protein", "Number of annotations", "Criteria")
#set criteria number for each
for(i in seq(1, nrow(CC_data_frame)-2, 3)) {
  CC_data_frame$Criteria[i] = 1
  CC_data_frame$Criteria[i+1] = 2
  CC_data_frame$Criteria[i+2] = 3
}
# set protein name for each
Protein_names = c("TNFAIP3", "BCL10", "N4BP1", "RC3H1", "RC3H2", "MALT1", "LIMA1", "MAP3K14", "RELB", "RBCK1", "CYLD", "ZC3H12A")
for(i in seq(1, nrow(CC_data_frame), 3)) {
  constant = (i + 2) / 3  
  CC_data_frame$Protein[i] = Protein_names[constant]
  CC_data_frame$Protein[i+1] = Protein_names[constant]
  CC_data_frame$Protein[i+2] = Protein_names[constant]
}
# now addnumber of annotations calculated previously
annotations = c(TNFAIP3_criteria_1, TNFAIP3_criteria_2, TNFAIP3_criteria_3, 
                BCL10_criteria_1, BCL10_criteria_2, BCL10_criteria_3, 
                N4BP1_criteria_1, N4BP1_criteria_2, N4BP1_criteria_3, 
                RC3H1_criteria_1, RC3H1_criteria_2, RC3H1_criteria_3,
                RC3H2_criteria_1, RC3H2_criteria_2, RC3H2_criteria_3, 
                MALT1_criteria_1, MALT1_criteria_2, MALT1_criteria_3, 
                LIMA1_criteria_1, LIMA1_criteria_2, LIMA1_criteria_3, 
                MAP3K14_criteria_1, MAP3K14_criteria_2, MAP3K14_criteria_3,
                RELB_criteria_1, RELB_criteria_2, RELB_criteria_3, 
                RBCK1_criteria_1, RBCK1_criteria_2, RBCK1_criteria_3, 
                CYLD_criteria_1, CYLD_criteria_2, CYLD_criteria_3, 
                ZC3H12A_criteria_1, ZC3H12A_criteria_2, ZC3H12A_criteria_3)
# add info to dataframe
CC_data_frame[2] = annotations 
colnames(CC_data_frame)[2] = "Number_of_annotations"
CC_data_frame[2] = as.numeric(unlist(CC_data_frame[2]))
CC_data_frame[3] = as.character(unlist(CC_data_frame[3]))
# now plot with ggplot2
library(ggplot2)
ggplot(CC_data_frame, aes(x=Protein, y=Number_of_annotations, fill=Criteria)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_y_continuous(limits = c(0,3), breaks = c(0, 1, 2, 3), labels = c(0, 1, 2, 3)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x=element_text(size = 12, face = "bold")) +
  theme(axis.title.y=element_text(size = 12, face = "bold")) +
  theme(axis.text.x=element_text(size = 12, face = "bold")) +
  theme(axis.text.y=element_text(size = 12, face = "bold")) +
  theme(axis.line.y = element_line(color = "black",size = 0.5)) +
  theme(legend.text = element_text(size = 12, face = "bold")) +
  theme(legend.title = element_text(size = 12, face = "bold")) + 
  theme(axis.ticks.x = element_blank())
ggsave(width = 30, height = 7, units = "cm", filename ="Supplementary_Fig2_GO_CC.pdf")



