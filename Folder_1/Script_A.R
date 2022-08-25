##### SCRIPT A #####

##### R script to compute GO-2-Substrates rank of proteins in PSSM-winnowing mode #####

# load required packages
library("dplyr")
library("stringr")

# define not in function
'%!in%' <- function(x,y)!('%in%'(x,y))

# import results from uniprotKB
output_uniprotKB = read.delim(file="uniprot-proteome_UP000005640+reviewed_yes.tab", header=TRUE, sep = "\t")

# import cutsites
cutsites = read.csv(file="fimo_sensitivity_1.csv", header=TRUE)
cutsites = data.frame(cutsites)

# rename columns in cutsites: rename sequence_name to Entry.name to avoid confusion
names(cutsites)[names(cutsites) == 'sequence_name'] = 'Entry.name'
names(cutsites)[names(cutsites) == 'matched_sequence'] = 'Matched.sequence'

# for convenience / compatiblity with previous code, make new frame called sequences
sequences = output_uniprotKB

# import known substrates
known_substrates = read.csv(file="known_substrates.csv", header=FALSE)
known_substrates_entry_name = as.vector(known_substrates[,2]) 
known_substrates = as.vector(known_substrates[,1])

# subset known substrates from output_uniprotKB
output_uniprotKB_known = output_uniprotKB[output_uniprotKB$Entry %in% known_substrates,]

# just in case there is something in work environment called 'c', remove c for do.call function that follows
remove(c)
# obtain list of unique biological process annotations for all MALT1 known substrates
output_uniprotKB_known_BP = str_extract_all(output_uniprotKB_known$Gene.ontology..biological.process., "\\[\\S*\\]")
output_uniprotKB_known_BP = do.call(c, output_uniprotKB_known_BP)
output_uniprotKB_known_BP = str_remove(output_uniprotKB_known_BP, "\\[")
output_uniprotKB_known_BP = str_remove(output_uniprotKB_known_BP, "\\]")

# obtain list & output_uniprotKB indices of unique biological process annotations 
# enriched in MALT1 known substrates, in different 'shells' (criteria) 
new_GO_BP_string = read.csv(file="GO_BP_string.csv", header=TRUE)
ordered_GO_BP_string = new_GO_BP_string[order(new_GO_BP_string$strength, decreasing = TRUE),]
ordered_GO_BP_string$strength_transform = 10^(ordered_GO_BP_string$strength)

# starting with last 5% rows, iteritively subset last n number of rows
# work out linear model for these data points, and work out Rsquared
# put in a data frame
k = nrow(ordered_GO_BP_string)

GO_BP_linear_models = data.frame(matrix(ncol=2,nrow=k))
GO_BP_linear_models$gene_count = NA
GO_BP_linear_models$rsquared_adj = NA
GO_BP_iteration = ls()

temp_index = 1
temp_constant = round(k*0.05)
# use worst 5% values as 'base' of line, then skip 5% of values,
# then iteratively add values, and see which has best line of best fit
# i.e attempt to find best linear model to describe the points that are least informative
for(i in 1:k-temp_constant) {
  GO_BP_iteration = ordered_GO_BP_string[(k-temp_constant-i):k,]
  strength_temp = GO_BP_iteration$strength
  gene_count_temp = GO_BP_iteration$background.gene.count
  summary_temp =  summary(lm(strength_temp ~ gene_count_temp))
  rsquared_temp = summary_temp$adj.r.squared
  intercept_temp = summary_temp$coefficients[1,1]
  gradient_temp = summary_temp$coefficients[2,1]
  GO_BP_linear_models$gene_count[temp_index] = min(GO_BP_iteration$background.gene.count)
  GO_BP_linear_models$gene_count_index[temp_index] = as.numeric(k-temp_constant-i)
  GO_BP_linear_models$rsquared_adj[temp_index] = rsquared_temp
  GO_BP_linear_models$intercept[temp_index] = intercept_temp
  GO_BP_linear_models$gradient[temp_index] = gradient_temp
  temp_index = temp_index + 1
}

# now trim away top 5% of rows (as they will likely be among best R^2)
GO_BP_linear_models = GO_BP_linear_models[-1:-temp_constant,]
# obtain intercept at max_rsquared_adj from GO_BP_linear_models
intercept_GO_BP = GO_BP_linear_models$intercept[which.max(GO_BP_linear_models$rsquared_adj)]
# obtain average 
intercept_GO_BP_filtered = ordered_GO_BP_string %>% filter(strength>=intercept_GO_BP)
# define list to be used as shell3 for GO_BP
GO_BP_shell3 = ordered_GO_BP_string %>% filter(strength<intercept_GO_BP)
# define list to be used as shell1 for GO_BP
temp_constant = median(intercept_GO_BP_filtered$strength)
GO_BP_shell1 = intercept_GO_BP_filtered %>% filter(strength>=temp_constant)
# define list to be used as shell2 for GO_BP
GO_BP_shell2 = intercept_GO_BP_filtered
GO_BP_shell1 = as.vector(GO_BP_shell1[,1])
GO_BP_shell2 = as.vector(GO_BP_shell2[,1])
GO_BP_shell3 = as.vector(GO_BP_shell3[,1])

# add shell 1 info onto output_uniprotKB
output_uniprotKB_GO_BP = str_extract_all(output_uniprotKB$Gene.ontology..biological.process., "\\[\\S*\\]")
for(i in 1:length(output_uniprotKB_GO_BP)) {
  output_uniprotKB_GO_BP[[i]] = str_remove(output_uniprotKB_GO_BP[[i]], "\\[")
  output_uniprotKB_GO_BP[[i]] = str_remove(output_uniprotKB_GO_BP[[i]], "\\]")
}

output_uniprotKB_GO_BP_index = NA
for(i in 1:length(output_uniprotKB_GO_BP)) {
  k = output_uniprotKB_GO_BP[[i]] %in% GO_BP_shell1
  if("TRUE" %in% k) {
    output_uniprotKB_GO_BP_index[i] = TRUE
  }
  if(TRUE %!in% k) {
    output_uniprotKB_GO_BP_index[i] = FALSE
  }
  remove(k)
  print(i)
}
remove(k)
output_uniprotKB_GO_BP_shell1 = output_uniprotKB[which(output_uniprotKB_GO_BP_index == TRUE),] 
output_uniprotKB_GO_BP_shell1_indices = which(output_uniprotKB_GO_BP_index == TRUE)


# shell 2_GO_BP
output_uniprotKB_GO_BP_index = NA
for(i in 1:length(output_uniprotKB_GO_BP)) {
  k = output_uniprotKB_GO_BP[[i]] %in% GO_BP_shell2
  if("TRUE" %in% k) {
    output_uniprotKB_GO_BP_index[i] = TRUE
  }
  if(TRUE %!in% k) {
    output_uniprotKB_GO_BP_index[i] = FALSE
  }
  remove(k)
  print(i)
}
remove(k)
output_uniprotKB_GO_BP_shell2 = output_uniprotKB[which(output_uniprotKB_GO_BP_index == TRUE),] 
output_uniprotKB_GO_BP_shell2_indices = which(output_uniprotKB_GO_BP_index == TRUE)

# shell 3_GO_BP
output_uniprotKB_GO_BP_index = NA
for(i in 1:length(output_uniprotKB_GO_BP)) {
  k = output_uniprotKB_GO_BP[[i]] %in% GO_BP_shell3
  if("TRUE" %in% k) {
    output_uniprotKB_GO_BP_index[i] = TRUE
  }
  if(TRUE %!in% k) {
    output_uniprotKB_GO_BP_index[i] = FALSE
  }
  remove(k)
  print(i)
}
remove(k)
output_uniprotKB_GO_BP_shell3 = output_uniprotKB[which(output_uniprotKB_GO_BP_index == TRUE),] 
output_uniprotKB_GO_BP_shell3_indices = which(output_uniprotKB_GO_BP_index == TRUE)


# add in columns to output_uniprotKB, append GO_BP_shell info 
output_uniprotKB$GO_BP_shell1 = NA
output_uniprotKB$GO_BP_shell1[output_uniprotKB_GO_BP_shell1_indices] = 1
output_uniprotKB$GO_BP_shell2 = NA
output_uniprotKB$GO_BP_shell2[output_uniprotKB_GO_BP_shell2_indices] = 1
output_uniprotKB$GO_BP_shell3 = NA
output_uniprotKB$GO_BP_shell3[output_uniprotKB_GO_BP_shell3_indices] = 1


#### Repeat for Molecular Function
new_GO_MF_string = read.csv(file="GO_MF_string.csv", header=TRUE)
ordered_GO_MF_string = new_GO_MF_string[order(new_GO_MF_string$strength, decreasing = TRUE),]
ordered_GO_MF_string$strength_transform = 10^(ordered_GO_MF_string$strength)

# starting with last 5% rows, iteritively subset last n number of rows
# work out linear model for these data points, and work out Rsquared
# put in a data frame
k = nrow(ordered_GO_MF_string)

GO_MF_linear_models = data.frame(matrix(ncol=2,nrow=k))
GO_MF_linear_models$gene_count = NA
GO_MF_linear_models$rsquared_adj = NA
GO_MF_iteration = ls()

temp_index = 1
temp_constant = round(k*0.05)
# use worst 5% values as 'base' of line, then skip 5% of values,
# then iteratively add values, and see which has best line of best fit
# i.e attempt to find best linear model to describe the points that we don't want!
for(i in 1:k-temp_constant) {
  GO_MF_iteration = ordered_GO_MF_string[(k-temp_constant-i):k,]
  strength_temp = GO_MF_iteration$strength
  gene_count_temp = GO_MF_iteration$background.gene.count
  summary_temp =  summary(lm(strength_temp ~ gene_count_temp))
  rsquared_temp = summary_temp$adj.r.squared
  intercept_temp = summary_temp$coefficients[1,1]
  gradient_temp = summary_temp$coefficients[2,1]
  GO_MF_linear_models$gene_count[temp_index] = min(GO_MF_iteration$background.gene.count)
  GO_MF_linear_models$gene_count_index[temp_index] = as.numeric(k-temp_constant-i)
  GO_MF_linear_models$rsquared_adj[temp_index] = rsquared_temp
  GO_MF_linear_models$intercept[temp_index] = intercept_temp
  GO_MF_linear_models$gradient[temp_index] = gradient_temp
  temp_index = temp_index + 1
}

# now trim away top 5% of rows (as they will likely be among best R^2)
GO_MF_linear_models = GO_MF_linear_models[-1:-temp_constant,]
# obtain intercept at max_rsquared_adj from GO_MF_linear_models
intercept_GO_MF = GO_MF_linear_models$intercept[which.max(GO_MF_linear_models$rsquared_adj)]
# obtain average 
intercept_GO_MF_filtered = ordered_GO_MF_string %>% filter(strength>=intercept_GO_MF)
# define list to be used as shell3 for GO_MF
GO_MF_shell3 = ordered_GO_MF_string %>% filter(strength<intercept_GO_MF)
# define list to be used as shell1 for GO_MF
temp_constant = median(intercept_GO_MF_filtered$strength)
GO_MF_shell1 = intercept_GO_MF_filtered %>% filter(strength>=temp_constant)
# define list to be used as shell1 for GO_MF
GO_MF_shell2 = intercept_GO_MF_filtered
GO_MF_shell1 = as.vector(GO_MF_shell1[,1])
GO_MF_shell2 = as.vector(GO_MF_shell2[,1])
GO_MF_shell3 = as.vector(GO_MF_shell3[,1])

# add shell 1 info onto output_uniprotKB
output_uniprotKB_GO_MF = str_extract_all(output_uniprotKB$Gene.ontology..molecular.function., "\\[\\S*\\]")
for(i in 1:length(output_uniprotKB_GO_MF)) {
  output_uniprotKB_GO_MF[[i]] = str_remove(output_uniprotKB_GO_MF[[i]], "\\[")
  output_uniprotKB_GO_MF[[i]] = str_remove(output_uniprotKB_GO_MF[[i]], "\\]")
}

output_uniprotKB_GO_MF_index = NA
for(i in 1:length(output_uniprotKB_GO_MF)) {
  k = output_uniprotKB_GO_MF[[i]] %in% GO_MF_shell1
  if("TRUE" %in% k) {
    output_uniprotKB_GO_MF_index[i] = TRUE
  }
  if(TRUE %!in% k) {
    output_uniprotKB_GO_MF_index[i] = FALSE
  }
  remove(k)
  print(i)
}
remove(k)
output_uniprotKB_GO_MF_shell1 = output_uniprotKB[which(output_uniprotKB_GO_MF_index == TRUE),] 
output_uniprotKB_GO_MF_shell1_indices = which(output_uniprotKB_GO_MF_index == TRUE)


# shell 2_GO_MF
output_uniprotKB_GO_MF_index = NA
for(i in 1:length(output_uniprotKB_GO_MF)) {
  k = output_uniprotKB_GO_MF[[i]] %in% GO_MF_shell2
  if("TRUE" %in% k) {
    output_uniprotKB_GO_MF_index[i] = TRUE
  }
  if(TRUE %!in% k) {
    output_uniprotKB_GO_MF_index[i] = FALSE
  }
  remove(k)
  print(i)
}
remove(k)
output_uniprotKB_GO_MF_shell2 = output_uniprotKB[which(output_uniprotKB_GO_MF_index == TRUE),] 
output_uniprotKB_GO_MF_shell2_indices = which(output_uniprotKB_GO_MF_index == TRUE)

# shell 3_GO_MF
output_uniprotKB_GO_MF_index = NA
for(i in 1:length(output_uniprotKB_GO_MF)) {
  k = output_uniprotKB_GO_MF[[i]] %in% GO_MF_shell3
  if("TRUE" %in% k) {
    output_uniprotKB_GO_MF_index[i] = TRUE
  }
  if(TRUE %!in% k) {
    output_uniprotKB_GO_MF_index[i] = FALSE
  }
  remove(k)
  print(i)
}
remove(k)
output_uniprotKB_GO_MF_shell3 = output_uniprotKB[which(output_uniprotKB_GO_MF_index == TRUE),] 
output_uniprotKB_GO_MF_shell3_indices = which(output_uniprotKB_GO_MF_index == TRUE)


# add in columns to output_uniprotKB, append GO_MF_shell info 
output_uniprotKB$GO_MF_shell1 = NA
output_uniprotKB$GO_MF_shell1[output_uniprotKB_GO_MF_shell1_indices] = 1
output_uniprotKB$GO_MF_shell2 = NA
output_uniprotKB$GO_MF_shell2[output_uniprotKB_GO_MF_shell2_indices] = 1
output_uniprotKB$GO_MF_shell3 = NA
output_uniprotKB$GO_MF_shell3[output_uniprotKB_GO_MF_shell3_indices] = 1


#### Repeat for Cellular Compartment
new_GO_CC_string = read.csv(file="GO_CC_string.csv", header=TRUE)
ordered_GO_CC_string = new_GO_CC_string[order(new_GO_CC_string$strength, decreasing = TRUE),]
ordered_GO_CC_string$strength_transform = 10^(ordered_GO_CC_string$strength)

# starting with last 5% rows, iteritively subset last n number of rows
# work out linear model for these data points, and work out Rsquared
# put in a data frame
k = nrow(ordered_GO_CC_string)

GO_CC_linear_models = data.frame(matrix(ncol=2,nrow=k))
GO_CC_linear_models$gene_count = NA
GO_CC_linear_models$rsquared_adj = NA
GO_CC_iteration = ls()

temp_index = 1
temp_constant = round(k*0.05)
# use worst 5% values as 'base' of line, then skip 5% of values,
# then iteratively add values, and see which has best line of best fit
# i.e attempt to find best linear model to describe the points that we don't want!
for(i in 1:k-temp_constant) {
  GO_CC_iteration = ordered_GO_CC_string[(k-temp_constant-i):k,]
  strength_temp = GO_CC_iteration$strength
  gene_count_temp = GO_CC_iteration$background.gene.count
  summary_temp =  summary(lm(strength_temp ~ gene_count_temp))
  rsquared_temp = summary_temp$adj.r.squared
  intercept_temp = summary_temp$coefficients[1,1]
  gradient_temp = summary_temp$coefficients[2,1]
  GO_CC_linear_models$gene_count[temp_index] = min(GO_CC_iteration$background.gene.count)
  GO_CC_linear_models$gene_count_index[temp_index] = as.numeric(k-temp_constant-i)
  GO_CC_linear_models$rsquared_adj[temp_index] = rsquared_temp
  GO_CC_linear_models$intercept[temp_index] = intercept_temp
  GO_CC_linear_models$gradient[temp_index] = gradient_temp
  temp_index = temp_index + 1
}

# now trim away top 5% of rows (as they will likely be among best R^2)
GO_CC_linear_models = GO_CC_linear_models[-1:-temp_constant,]
# obtain intercept at max_rsquared_adj from GO_CC_linear_models
intercept_GO_CC = GO_CC_linear_models$intercept[which.max(GO_CC_linear_models$rsquared_adj)]
# obtain average 
intercept_GO_CC_filtered = ordered_GO_CC_string %>% filter(strength>=intercept_GO_CC)
# define list to be used as shell3 for GO_CC
GO_CC_shell3 = ordered_GO_CC_string %>% filter(strength<intercept_GO_CC)
# define list to be used as shell1 for GO_CC
temp_constant = median(intercept_GO_CC_filtered$strength)
GO_CC_shell1 = intercept_GO_CC_filtered %>% filter(strength>=temp_constant)
# define list to be used as shell1 for GO_CC
GO_CC_shell2 = intercept_GO_CC_filtered
GO_CC_shell1 = as.vector(GO_CC_shell1[,1])
GO_CC_shell2 = as.vector(GO_CC_shell2[,1])
GO_CC_shell3 = as.vector(GO_CC_shell3[,1])


# add shell 1 info onto output_uniprotKB
output_uniprotKB_GO_CC = str_extract_all(output_uniprotKB$Gene.ontology..cellular.component., "\\[\\S*\\]")
for(i in 1:length(output_uniprotKB_GO_CC)) {
  output_uniprotKB_GO_CC[[i]] = str_remove(output_uniprotKB_GO_CC[[i]], "\\[")
  output_uniprotKB_GO_CC[[i]] = str_remove(output_uniprotKB_GO_CC[[i]], "\\]")
}

output_uniprotKB_GO_CC_index = NA
for(i in 1:length(output_uniprotKB_GO_CC)) {
  k = output_uniprotKB_GO_CC[[i]] %in% GO_CC_shell1
  if("TRUE" %in% k) {
    output_uniprotKB_GO_CC_index[i] = TRUE
  }
  if(TRUE %!in% k) {
    output_uniprotKB_GO_CC_index[i] = FALSE
  }
  remove(k)
  print(i)
}
remove(k)
output_uniprotKB_GO_CC_shell1 = output_uniprotKB[which(output_uniprotKB_GO_CC_index == TRUE),] 
output_uniprotKB_GO_CC_shell1_indices = which(output_uniprotKB_GO_CC_index == TRUE)


# shell 2_GO_CC
output_uniprotKB_GO_CC_index = NA
for(i in 1:length(output_uniprotKB_GO_CC)) {
  k = output_uniprotKB_GO_CC[[i]] %in% GO_CC_shell2
  if("TRUE" %in% k) {
    output_uniprotKB_GO_CC_index[i] = TRUE
  }
  if(TRUE %!in% k) {
    output_uniprotKB_GO_CC_index[i] = FALSE
  }
  remove(k)
  print(i)
}
remove(k)
output_uniprotKB_GO_CC_shell2 = output_uniprotKB[which(output_uniprotKB_GO_CC_index == TRUE),] 
output_uniprotKB_GO_CC_shell2_indices = which(output_uniprotKB_GO_CC_index == TRUE)

# shell 3_GO_CC
output_uniprotKB_GO_CC_index = NA
for(i in 1:length(output_uniprotKB_GO_CC)) {
  k = output_uniprotKB_GO_CC[[i]] %in% GO_CC_shell3
  if("TRUE" %in% k) {
    output_uniprotKB_GO_CC_index[i] = TRUE
  }
  if(TRUE %!in% k) {
    output_uniprotKB_GO_CC_index[i] = FALSE
  }
  remove(k)
  print(i)
}
remove(k)
output_uniprotKB_GO_CC_shell3 = output_uniprotKB[which(output_uniprotKB_GO_CC_index == TRUE),] 
output_uniprotKB_GO_CC_shell3_indices = which(output_uniprotKB_GO_CC_index == TRUE)

# add in columns to output_uniprotKB, append GO_CC_shell info 
output_uniprotKB$GO_CC_shell1 = NA
output_uniprotKB$GO_CC_shell1[output_uniprotKB_GO_CC_shell1_indices] = 1
output_uniprotKB$GO_CC_shell2 = NA
output_uniprotKB$GO_CC_shell2[output_uniprotKB_GO_CC_shell2_indices] = 1
output_uniprotKB$GO_CC_shell3 = NA
output_uniprotKB$GO_CC_shell3[output_uniprotKB_GO_CC_shell3_indices] = 1

#### InterPro 
# subset rows of output_uniprotKB containing interpro_shell1 terms
interpro_shell1 = read.csv(file="InterPro_shell1.csv", header=TRUE)
interpro_shell1 = as.vector(interpro_shell1[,1])
output_uniprotKB_interpro = str_split(output_uniprotKB$Cross.reference..InterPro, ";")
# to do this, obtain indices of matching terms, then subset
output_uniprotKB_interpro_index = NA
for(i in 1:length(output_uniprotKB_interpro)) {
  k = output_uniprotKB_interpro[[i]] %in% interpro_shell1
  if("TRUE" %in% k) {
    output_uniprotKB_interpro_index[i] = TRUE
  }
  if(TRUE %!in% k) {
    output_uniprotKB_interpro_index[i] = FALSE
  }
  remove(k)
  print(i)
}
remove(k)
output_uniprotKB_interpro_shell1 = output_uniprotKB[which(output_uniprotKB_interpro_index == TRUE),] 
output_uniprotKB_interpro_shell1_indices = which(output_uniprotKB_interpro_index == TRUE)

## Derive interpro_shell2 terms

# first extract known MALT1 substrates
genes_of_interest = read.csv("known_substrates.csv", header = FALSE)

# extract MALT1 substrates from output_uniprotKB
temp_interpro_table = output_uniprotKB[which(output_uniprotKB$Entry.name %in% genes_of_interest$V2),]
# extract all unique interpro terms for MALT1 substrates
temp_interpro = paste(temp_interpro_table$Cross.reference..InterPro., sep = ";", collapse = "")
temp_interpro = str_split(temp_interpro, ";", n = Inf, simplify = FALSE)
temp_interpro_unique = unique(temp_interpro)
temp_interpro_unique_column = matrix(temp_interpro_unique[[1]])
temp_interpro_unique_column = as.data.frame(temp_interpro_unique_column)
temp_interpro_unique_column = temp_interpro_unique_column %>% filter(V1 != "")
write.csv(temp_interpro_unique_column, "interpro_shell2.csv")
interpro_shell2 = temp_interpro_unique_column

# subset rows of output_uniprotKB containing interpro_shell2 terms
interpro_shell2 = as.vector(interpro_shell2$V1)
output_uniprotKB_interpro = str_split(output_uniprotKB$Cross.reference..InterPro, ";")
# to do this, obtain indices ofmatching terms, then subset
output_uniprotKB_interpro_index = NA
for(i in 1:length(output_uniprotKB_interpro)) {
  k = output_uniprotKB_interpro[[i]] %in% interpro_shell2
  if("TRUE" %in% k) {
    output_uniprotKB_interpro_index[i] = TRUE
  }
  if(TRUE %!in% k) {
    output_uniprotKB_interpro_index[i] = FALSE
  }
  remove(k)
  print(i)
}
output_uniprotKB_interpro_shell2 = output_uniprotKB[which(output_uniprotKB_interpro_index == TRUE),] 
output_uniprotKB_interpro_shell2_indices = which(output_uniprotKB_interpro_index == TRUE)

# add in columns to output_uniprotKB, append interpro_shell info 
output_uniprotKB$interpro_shell1 = NA
output_uniprotKB$interpro_shell1[output_uniprotKB_interpro_shell1_indices] = 1
output_uniprotKB$interpro_shell2 = NA
output_uniprotKB$interpro_shell2[output_uniprotKB_interpro_shell2_indices] = 1

##### PPIs
# import known interactors from BioGrid
biogrid = read.csv(file="BIOGRID-ORGANISM-Homo_sapiens-4.2.191.csv", header=TRUE)
# obtain list of proteins known to interact with MALT1
biogrid_MALT1 = biogrid[biogrid$Official.Symbol.Interactor.A %in% "MALT1",]
interacts_shell1 = unique(biogrid_MALT1$Official.Symbol.Interactor.B)
# obtain list of proteins known to interact with substrates of malt1; include gene name synonyms
# just in case there is something in work environment called 'c', remove c for do.call function that follows
remove(c)
known_gene_names = output_uniprotKB_known$Gene.names
known_gene_names = str_split(known_gene_names, " ")
known_gene_names = do.call(c, known_gene_names)
biogrid_substrates = biogrid[biogrid$Official.Symbol.Interactor.A %in% known_gene_names,]
interacts_shell2 = unique(biogrid_substrates$Official.Symbol.Interactor.B)


# subset rows of output_uniprotKB containing interacts_shell1 terms
output_uniprotKB_gene_names = str_split(output_uniprotKB$Gene.names, " ")
# to do this, obtain indeces ofmatching terms, then subset
output_uniprotKB_interacts_index = NA
for(i in 1:length(output_uniprotKB_gene_names)) {
  k = output_uniprotKB_gene_names[[i]] %in% interacts_shell1
  if("TRUE" %in% k) {
    output_uniprotKB_interacts_index[i] = TRUE
  }
  if(TRUE %!in% k) {
    output_uniprotKB_interacts_index[i] = FALSE
  }
  remove(k)
  print(i)
}
remove(k)
output_uniprotKB_interacts_shell1 = output_uniprotKB[which(output_uniprotKB_interacts_index == TRUE),] 
output_uniprotKB_interacts_shell1_indices = which(output_uniprotKB_interacts_index == TRUE)

# subset rows of output_uniprotKB containing interacts_shell2 terms
output_uniprotKB_gene_names = str_split(output_uniprotKB$Gene.names, " ")
# to do this, obtain indeces ofmatching terms, then subset
output_uniprotKB_interacts_index = NA
for(i in 1:length(output_uniprotKB_gene_names)) {
  k = output_uniprotKB_gene_names[[i]] %in% interacts_shell2
  if("TRUE" %in% k) {
    output_uniprotKB_interacts_index[i] = TRUE
  }
  if(TRUE %!in% k) {
    output_uniprotKB_interacts_index[i] = FALSE
  }
  remove(k)
  print(i)
}
remove(k)
output_uniprotKB_interacts_shell2 = output_uniprotKB[which(output_uniprotKB_interacts_index == TRUE),] 
output_uniprotKB_interacts_shell2_indices = which(output_uniprotKB_interacts_index == TRUE)

# add in columns to output_uniprotKB, append interacts_shell info 
output_uniprotKB$interacts_shell1 = NA
output_uniprotKB$interacts_shell1[output_uniprotKB_interacts_shell1_indices] = 1
output_uniprotKB$interacts_shell2 = NA
output_uniprotKB$interacts_shell2[output_uniprotKB_interacts_shell2_indices] = 1


#### Calculate scores for each feature and criteria

# first determine how many instances in proteome
GO_BP_shell1_proteome=length(which(output_uniprotKB$GO_BP_shell1==TRUE))
GO_BP_shell2_proteome=length(which(output_uniprotKB$GO_BP_shell2==TRUE))
GO_BP_shell3_proteome=length(which(output_uniprotKB$GO_BP_shell3==TRUE))
GO_MF_shell1_proteome=length(which(output_uniprotKB$GO_MF_shell1==TRUE))
GO_MF_shell2_proteome=length(which(output_uniprotKB$GO_MF_shell2==TRUE))
GO_MF_shell3_proteome=length(which(output_uniprotKB$GO_MF_shell3==TRUE))
GO_CC_shell1_proteome=length(which(output_uniprotKB$GO_CC_shell1==TRUE))
GO_CC_shell2_proteome=length(which(output_uniprotKB$GO_CC_shell2==TRUE))
GO_CC_shell3_proteome=length(which(output_uniprotKB$GO_CC_shell3==TRUE))
interpro_shell1_proteome=length(which(output_uniprotKB$interpro_shell1==TRUE))
interpro_shell2_proteome= length(which(output_uniprotKB$interpro_shell2==TRUE))
interacts_shell1_proteome=length(which(output_uniprotKB$interacts_shell1==TRUE))
interacts_shell2_proteome=length(which(output_uniprotKB$interacts_shell2==TRUE))

# next determine how many known MALT1 substrates meet this criteria
output_uniprotKB_known_substrates = output_uniprotKB[which(output_uniprotKB$Entry %in% known_substrates),]
number_known_substrates = nrow(output_uniprotKB_known_substrates)

BP_1_constant = length(which(output_uniprotKB_known_substrates$GO_BP_shell1 == 1)) / number_known_substrates
BP_2_constant = length(which(output_uniprotKB_known_substrates$GO_BP_shell2 == 1)) / number_known_substrates
BP_3_constant = length(which(output_uniprotKB_known_substrates$GO_BP_shell3 == 1)) / number_known_substrates
MF_1_constant = length(which(output_uniprotKB_known_substrates$GO_MF_shell1 == 1)) / number_known_substrates
MF_2_constant = length(which(output_uniprotKB_known_substrates$GO_MF_shell2 == 1)) / number_known_substrates
MF_3_constant = length(which(output_uniprotKB_known_substrates$GO_MF_shell3 == 1)) / number_known_substrates
CC_1_constant = length(which(output_uniprotKB_known_substrates$GO_CC_shell1 == 1)) / number_known_substrates
CC_2_constant = length(which(output_uniprotKB_known_substrates$GO_CC_shell2 == 1)) / number_known_substrates
CC_3_constant = length(which(output_uniprotKB_known_substrates$GO_CC_shell3 == 1)) / number_known_substrates
interpro_1_constant = length(which(output_uniprotKB_known_substrates$interpro_shell1 == 1)) / number_known_substrates
interpro_2_constant = length(which(output_uniprotKB_known_substrates$interpro_shell2 == 1)) / number_known_substrates
interacts_1_constant = length(which(output_uniprotKB_known_substrates$interacts_shell1 == 1)) / number_known_substrates
interacts_2_constant = length(which(output_uniprotKB_known_substrates$interacts_shell2 == 1)) / number_known_substrates


# next express this as a fraction and use to
# finally, normalise scores
output_uniprotKB$GO_BP_shell1_score= (output_uniprotKB$GO_BP_shell1 / GO_BP_shell1_proteome) * BP_1_constant 
output_uniprotKB$GO_BP_shell2_score= (output_uniprotKB$GO_BP_shell2 / GO_BP_shell2_proteome) * BP_2_constant 
output_uniprotKB$GO_BP_shell3_score= (output_uniprotKB$GO_BP_shell3 / GO_BP_shell3_proteome) * BP_3_constant 
output_uniprotKB$GO_MF_shell1_score= (output_uniprotKB$GO_MF_shell1 / GO_MF_shell1_proteome) * MF_1_constant 
output_uniprotKB$GO_MF_shell2_score= (output_uniprotKB$GO_MF_shell2 / GO_MF_shell2_proteome) * MF_2_constant 
output_uniprotKB$GO_MF_shell3_score= (output_uniprotKB$GO_MF_shell3 / GO_MF_shell3_proteome) * MF_3_constant 
output_uniprotKB$GO_CC_shell1_score= (output_uniprotKB$GO_CC_shell1 / GO_CC_shell1_proteome) * CC_1_constant
output_uniprotKB$GO_CC_shell2_score= (output_uniprotKB$GO_CC_shell2 / GO_CC_shell2_proteome) * CC_2_constant
output_uniprotKB$GO_CC_shell3_score= (output_uniprotKB$GO_CC_shell3 / GO_CC_shell3_proteome) * CC_3_constant
output_uniprotKB$interpro_shell1_score= (output_uniprotKB$interpro_shell1 / interpro_shell1_proteome) * interpro_1_constant 
output_uniprotKB$interpro_shell2_score= (output_uniprotKB$interpro_shell2 / interpro_shell2_proteome) * interpro_2_constant 
output_uniprotKB$interacts_shell1_score= (output_uniprotKB$interacts_shell1 / interacts_shell1_proteome) * interacts_1_constant
output_uniprotKB$interacts_shell2_score= (output_uniprotKB$interacts_shell2 / interacts_shell2_proteome) * interacts_2_constant

##### calculate combined Function score
for(i in 1:nrow(output_uniprotKB)) {
  output_uniprotKB$biological_combined_score[i] = sum(output_uniprotKB$GO_BP_shell1_score[i], output_uniprotKB$GO_BP_shell2_score[i], 
                                                      output_uniprotKB$GO_BP_shell3_score[i],
                                                      output_uniprotKB$GO_MF_shell1_score[i], output_uniprotKB$GO_MF_shell2_score[i], 
                                                      output_uniprotKB$GO_MF_shell3_score[i],
                                                      output_uniprotKB$GO_CC_shell1_score[i], output_uniprotKB$GO_CC_shell2_score[i],
                                                      output_uniprotKB$GO_CC_shell3_score[i],
                                                      output_uniprotKB$interpro_shell1_score[i], output_uniprotKB$interpro_shell2_score[i], output_uniprotKB$interacts_shell1_score[i],
                                                      output_uniprotKB$interacts_shell2_score[i], na.rm=TRUE)
}


# add empty column to cutsites to accomodate Sequences
cutsites$Gene_name=NA
cutsites$GO_BP_shell1_score=NA
cutsites$GO_BP_shell2_score=NA
cutsites$GO_BP_shell3_score=NA
cutsites$GO_MF_shell1_score=NA
cutsites$GO_MF_shell2_score=NA
cutsites$GO_MF_shell3_score=NA
cutsites$GO_CC_shell1_score=NA
cutsites$GO_CC_shell2_score=NA
cutsites$GO_CC_shell3_score=NA
cutsites$interpro_shell1_score=NA
cutsites$interpro_shell2_score=NA
cutsites$interacts_shell1_score=NA
cutsites$interacts_shell2_score=NA
cutsites$Protein_name = NA
cutsites$Sequence=NA
cutsites$biological_combined_score=NA


# loop, matching sequences$cutsites to sequences$Sequence
output_uniprotKB_in_cutsites = which(output_uniprotKB$Entry.name %in% cutsites$Entry.name)

for (i in 1:nrow(cutsites)) {
  for (j in output_uniprotKB_in_cutsites) {
    if(identical(as.vector(cutsites$Entry.name[i]), as.vector(output_uniprotKB$Entry.name[j]))) {
      a = as.vector(output_uniprotKB$Gene.names[j])
      b = as.vector(output_uniprotKB$GO_BP_shell1_score[j])
      c = as.vector(output_uniprotKB$GO_BP_shell2_score[j])
      d = as.vector(output_uniprotKB$GO_MF_shell1_score[j])
      e = as.vector(output_uniprotKB$GO_MF_shell2_score[j])
      f = as.vector(output_uniprotKB$GO_CC_shell1_score[j])
      g = as.vector(output_uniprotKB$GO_CC_shell2_score[j])
      h = as.vector(output_uniprotKB$interpro_shell1_score[j])
      k = as.vector(output_uniprotKB$interpro_shell2_score[j])
      l = as.vector(output_uniprotKB$interacts_shell1_score[j])
      m = as.vector(output_uniprotKB$interacts_shell2_score[j])
      n = as.vector(output_uniprotKB$Protein.names[j])
      o = as.vector(output_uniprotKB$Sequence[j])
      q = as.vector(output_uniprotKB$GO_BP_shell3_score[j])
      r = as.vector(output_uniprotKB$GO_MF_shell3_score[j])
      s = as.vector(output_uniprotKB$GO_CC_shell3_score[j])
      t = as.vector(output_uniprotKB$biological_combined_score[j])
      if(!is.na(a)) break }
  }
  cutsites$Gene_name[i]=a
  cutsites$GO_BP_shell1_score[i]=b
  cutsites$GO_BP_shell2_score[i]=c
  cutsites$GO_BP_shell3_score[i]=q
  cutsites$GO_MF_shell1_score[i]=d
  cutsites$GO_MF_shell2_score[i]=e
  cutsites$GO_MF_shell3_score[i]=r
  cutsites$GO_CC_shell1_score[i]=f
  cutsites$GO_CC_shell2_score[i]=g
  cutsites$GO_CC_shell3_score[i]=s
  cutsites$interpro_shell1_score[i]=h
  cutsites$interpro_shell2_score[i]= k
  cutsites$interacts_shell1_score[i]=l
  cutsites$interacts_shell2_score[i]=m
  cutsites$Protein_name[i] = n
  cutsites$Sequence[i]=o
  cutsites$biological_combined_score[i]=t
  a = NA; b=NA; c=NA; d=NA; e=NA; f=NA; g=NA; h=NA
  k=NA; l=NA; m=NA; n=NA; o=NA; q=NA; r=NA; s=NA; t=NA
  print(i)
}

# rank by biological_combined_score
cutsites$bio_rank = rank(-cutsites$biological_combined_score, ties.method = "min")

# order, so that 1st instance of each protein in cutsites is best bio_rank 
cutsites = cutsites[order(cutsites$bio_rank),]

# keep only best scoring bio score and make new data frame
cutsites_bio = cutsites[match(unique(cutsites$Entry.name), cutsites$Entry.name),]

# rank again by biological_combined_score
cutsites_bio$bio_rank = rank(-cutsites_bio$biological_combined_score, ties.method = "min")
a = NA
# now add cutsites_bio rank onto cutsites
for(i in 1:nrow(cutsites)) {
  for(j in 1:nrow(cutsites_bio)) {
    if(identical(as.vector(cutsites$Entry.name[i]), as.vector(cutsites_bio$Entry.name[j]))) {
      a = cutsites_bio$bio_rank[j]
      cutsites$bio_rank_unique[i] = a
      if(!is.na(a)) break 
    }
  }
  print(i)
}

cutsites$bio_rank = cutsites$bio_rank_unique



###### Sequence Module  #######

# load dplyr
library("dplyr")
# load stringr
library("stringr")

## cycle through output_uniprotKB, extracting rows where known substrate 
## shares P4-P1
cutsite_coop = read.csv(file="cleavage_sites.csv", header=TRUE)
cutsites$cutsite_coop = NA
#note - added in columns here to accomodate info later
cutsites$P6_P6prime_exact = NA
cutsites$P6_P2prime_exact = NA
cutsites$P6_P1prime_exact = NA
cutsites$P4_P2prime_exact = NA
cutsites$P4_P1prime_exact = NA


### when cutsite match is found, identify number of characters into $Sequences that match occurs
### create new column, containing all characters from N-1 to end of string

# add empty column to sequences to accomodate cutsite_coordinate (i.e. P1prime)
cutsites$cutsite_coordinate=NA

match_coordinate = c()
for (i in 1:nrow(cutsites)) {
  match_coordinate = str_locate(as.vector(cutsites$Sequence[i]), as.vector(cutsites$Matched.sequence[i]))
  # extract 'end' point of match sequence from sequence 
  match_coordinate = data.frame(match_coordinate)
  cutsite_coordinate = match_coordinate$end
  # append onto sequences data frame  
  cutsites$cutsite_coordinate[i] = cutsite_coordinate 
}

### extract sequence of putative neo-N-terminal peptide, assuming cleavage at cutsite
cutsites$neoN_term=NA
cutsites$neoN_term_length=NA
cutsites$neoN_term = substring(cutsites$Sequence, cutsites$cutsite_coordinate)
cutsites$neoN_term_length = (nchar(as.character(cutsites$Sequence)) - cutsites$cutsite_coordinate)+1

### now that sequences added to table, search for P6_P6prime exact etc
# to do this, need to generate sequences column for P6_P6 etc first to avoid pulling out
# all with the same Sequence
cutsites$P6_P6prime_sequence = substr(cutsites$Sequence, cutsites$cutsite_coordinate-6, cutsites$cutsite_coordinate+5)
cutsites$P6_P2prime_sequence = substr(cutsites$Sequence, cutsites$cutsite_coordinate-6, cutsites$cutsite_coordinate+1)
cutsites$P6_P1prime_sequence = substr(cutsites$Sequence, cutsites$cutsite_coordinate-6, cutsites$cutsite_coordinate)
cutsites$P4_P2prime_sequence = substr(cutsites$Sequence, cutsites$cutsite_coordinate-4, cutsites$cutsite_coordinate+1)
cutsites$P4_P1prime_sequence = substr(cutsites$Sequence, cutsites$cutsite_coordinate-4, cutsites$cutsite_coordinate)

cutsites$P6_P6prime_exact = NA
temp_frame = data.frame(matrix(, nrow=0, ncol=ncol(output_uniprotKB)))
colnames(temp_frame) = colnames(output_uniprotKB)
for (i in 1:nrow(cutsites)) {
  for (j in 1:nrow(cutsite_coop)) {
    if(grepl(as.vector(cutsite_coop$P6_P6prime[j]), as.vector(cutsites$P6_P6prime_sequence[i])))
      cutsites$P6_P6prime_exact[i] = 1
  }}
cutsite_P6_P6prime_only = cutsites %>% filter(P6_P6prime_exact == 1)

cutsites$P6_P2prime_exact = NA
temp_frame = data.frame(matrix(, nrow=0, ncol=ncol(output_uniprotKB)))
colnames(temp_frame) = colnames(output_uniprotKB)
for (i in 1:nrow(cutsites)) {
  for (j in 1:nrow(cutsite_coop)) {
    if(grepl(as.vector(cutsite_coop$P6_P2prime[j]), as.vector(cutsites$P6_P2prime_sequence[i])))
      cutsites$P6_P2prime_exact[i] = 1
  }}
cutsite_P6_P2prime_only = cutsites %>% filter(P6_P2prime_exact == 1)

cutsites$P6_P1prime_exact = NA
temp_frame = data.frame(matrix(, nrow=0, ncol=ncol(output_uniprotKB)))
colnames(temp_frame) = colnames(output_uniprotKB)
for (i in 1:nrow(cutsites)) {
  for (j in 1:nrow(cutsite_coop)) {
    if(grepl(as.vector(cutsite_coop$P6_P1prime[j]), as.vector(cutsites$P6_P1prime_sequence[i])))
      cutsites$P6_P1prime_exact[i] = 1
  }}
cutsite_P6_P1prime_only = cutsites %>% filter(P6_P1prime_exact == 1)

cutsites$P4_P2prime_exact = NA
temp_frame = data.frame(matrix(, nrow=0, ncol=ncol(output_uniprotKB)))
colnames(temp_frame) = colnames(output_uniprotKB)
for (i in 1:nrow(cutsites)) {
  for (j in 1:nrow(cutsite_coop)) {
    if(grepl(as.vector(cutsite_coop$P4_P2prime[j]), as.vector(cutsites$P4_P2prime_sequence[i])))
      cutsites$P4_P2prime_exact[i] = 1
  }}
cutsite_P4_P2prime_only = cutsites %>% filter(P4_P2prime_exact == 1)

cutsites$P4_P1prime_exact = NA
temp_frame = data.frame(matrix(, nrow=0, ncol=ncol(output_uniprotKB)))
colnames(temp_frame) = colnames(output_uniprotKB)
for (i in 1:nrow(cutsites)) {
  for (j in 1:nrow(cutsite_coop)) {
    if(grepl(as.vector(cutsite_coop$P4_P1prime[j]), as.vector(cutsites$P4_P1prime_sequence[i])))
      cutsites$P4_P1prime_exact[i] = 1
  }}
cutsite_P4_P1prime_only = cutsites %>% filter(P4_P1prime_exact == 1)
cutsite_P6_P6prime_only = cutsites %>% filter(P6_P6prime_exact == 1)
cutsite_P6_P2prime_only = cutsites %>% filter(P6_P2prime_exact == 1)
cutsite_P6_P1prime_only = cutsites %>% filter(P6_P1prime_exact == 1)
cutsite_P4_P2prime_only = cutsites %>% filter(P4_P2prime_exact == 1)


# add exact_match_scores
exact_match_shell1_indices = which(cutsites$P4_P2prime_exact == TRUE)
exact_match_shell2_indices = which(cutsites$P4_P1prime_exact == TRUE)

# add in columns to output_uniprotKB, append interacts_shell info 
cutsites$exact_match_shell1 = NA
cutsites$exact_match_shell1[exact_match_shell1_indices] = 1
cutsites$exact_match_shell2 = NA
cutsites$exact_match_shell2[exact_match_shell2_indices] = 1


# account for some proteins having multiple exact matches to known cutsites
# for each line in cutsites, extract all lines with same protein, and sum
# number of lines with exact_match_shell_1 or 2 == 1
exact_match_shell1_cutsites_subset = cutsites %>% filter(exact_match_shell1 == 1)
exact_match_shell2_cutsites_subset = cutsites %>% filter(exact_match_shell2 == 1)
cutsites$exact_match_shell1_sum = NA
cutsites$exact_match_shell2_sum = NA

# note - make sure that the specific cutsite being annotated is an exact match
# to avoid a non-exact match getting a score, just because an exact match
# exists elsewhere in protein
for(i in 1:nrow(cutsites)) {
  if(cutsites$Entry.name[i] %in% exact_match_shell1_cutsites_subset$Entry.name) {
    temp_frame = cutsites %>% filter(cutsites$Entry.name == cutsites$Entry.name[i])
    #updated here with loop
    if(cutsites$Matched.sequence[i] %in% cutsite_coop$P4_P1prime) {
      cutsites$exact_match_shell1_sum[i] = sum(temp_frame$exact_match_shell1, na.rm = TRUE)
    }
  }
}

for(i in 1:nrow(cutsites)) {
  if(cutsites$Entry.name[i] %in% exact_match_shell2_cutsites_subset$Entry.name) {
    temp_frame = cutsites %>% filter(cutsites$Entry.name == cutsites$Entry.name[i])
    cutsites$exact_match_shell2_sum[i] = sum(temp_frame$exact_match_shell2, na.rm = TRUE)
  }
}

### assign greater score for substrates with more exact matches to known cleavage site
cutsites$exact_match_shell1_score = NA
cutsites$exact_match_shell2_score = NA
for(i in 1:nrow(cutsites)) {
  if (!is.na(cutsites$exact_match_shell1[i])) {
    cutsites$exact_match_shell1_score[i] = cutsites$exact_match_shell1_sum[i] / length(which(cutsites$exact_match_shell1 > 0))
  }
  if (!is.na(cutsites$exact_match_shell2[i])) {
    cutsites$exact_match_shell2_score[i] = cutsites$exact_match_shell2_sum[i] / length(which(cutsites$exact_match_shell2 > 0))
  }
}

cutsites_backup = cutsites

######## ORTHOLOGS ###########
# run ortholog_script.R here
all_orthologs = read.csv("orthologs.csv", header = TRUE)

# subset 01:01 groups etc into separate tables
ortholog_01_01 = all_orthologs %>% filter(type == "01:01")
ortholog_n_01 = all_orthologs %>% filter(type == "n:1")
ortholog_01_m = all_orthologs %>% filter(type == "1:m")
ortholog_n_m = all_orthologs %>% filter(type == "n:m")

#NEW workflow!
# add on mouse entry_name, sequences and type to cutsites subtables
cutsites$Mouse_entry_name = NA
cutsites$Mouse_sequence = NA
cutsites$Ortholog_type = NA
temp_frame = data.frame(matrix(, nrow=0, ncol=ncol(cutsites)))
colnames(temp_frame) = colnames(cutsites)
# use the following index further on in code
all_orthologs_in_cutsites = which(all_orthologs$human_uniprot %in% cutsites$Entry.name)
for (j in 1:nrow(cutsites)) {
  for (i in all_orthologs_in_cutsites) {
    if(identical(as.vector(all_orthologs$human_uniprot[i]), as.character(cutsites$Entry.name[j]))) {
      cutsites$Mouse_entry_name[j] = as.character(all_orthologs$mouse_uniprot[i])
      cutsites$Mouse_sequence[j] = as.character(all_orthologs$mouse_sequence[i])
      cutsites$Ortholog_type[j] = as.character(all_orthologs$type[i])
      temp_frame = rbind(temp_frame, cutsites[j,]) 
      break }
  }
  print(j)
}
cutsites_all_orthologs = temp_frame
cutsites_all_orthologs = unique(cutsites_all_orthologs)

all_orthologs_human_ids = cutsites_all_orthologs$Entry.name   

# add column to show if mouse ortholog exists or not
cutsites$Mouse_ortholog_exists = NA
Mouse_ortholog_exists_index = which(!is.na(cutsites$Mouse_entry_name))
for(i in Mouse_ortholog_exists_index) {
  cutsites$Mouse_ortholog_exists[i] = 1
}

cutsites_no_ortholog = cutsites %>% filter(is.na(cutsites$Mouse_ortholog_exists))

# if mouse ortholog exists, test for conservation of cutsite human sequence in mouse sequence      
cutsites$Mouse_Human_exact_conserved = NA
for (i in 1:nrow(cutsites)) {
  if(grepl(as.vector(cutsites$Matched.sequence[i]), as.vector(cutsites$Mouse_sequence[i]))) {
    cutsites$Mouse_Human_exact_conserved[i] = "1"
  }
  print(i)
}

# test for exact match to a known malt1 substrate cutsite in mouse sequence
cutsites$Mouse_P4_P1prime_exact = NA
cutsites$Mouse_P4_P1prime_exact_sequence = NA
unique_cutsite_coop = unique(cutsite_coop$P4_P1prime)
for (i in 1:nrow(cutsites)) {
  for (j in 1:length(unique_cutsite_coop)) {
    if(grepl(as.vector(unique_cutsite_coop[j]), as.vector(cutsites$Mouse_sequence[i]))) {
      cutsites$Mouse_P4_P1prime_exact[i] = 1
      cutsites$Mouse_P4_P1prime_exact_sequence[i] = as.character(unique_cutsite_coop[j]) }   
  }
  print(i)
}

all_orthologs$P4_P1prime_exact_in_mouse = NA
all_orthologs$P4_P1prime_exact_in_human = NA
for(i in 1:nrow(all_orthologs)) {
  for(j in 1:length(unique_cutsite_coop)) {
    if(grepl(as.vector(unique_cutsite_coop[j]), as.vector(all_orthologs$mouse_sequence[i]))) {
      all_orthologs$P4_P1prime_exact_in_mouse[i] = 1
    }
    if(grepl(as.vector(unique_cutsite_coop[j]), as.vector(all_orthologs$human_sequence[i]))) {
      all_orthologs$P4_P1prime_exact_in_human[i] = 1
    }
  }
  print(i)
}

all_orthologs$both_species_P4_P1prime = NA
for(i in 1:nrow(all_orthologs)) {
  if(!is.na(all_orthologs$P4_P1prime_exact_in_human[i])) {
    if(all_orthologs$P4_P1prime_exact_in_human[i] == 1) {
      if(!is.na(all_orthologs$P4_P1prime_exact_in_mouse[i])) {
        if(all_orthologs$P4_P1prime_exact_in_mouse[i] == 1) {
          all_orthologs$both_species_P4_P1prime[i] = 1 } 
      } 
    }
  }
  print(i)
}

# append info onto cutsites
both_species_P4_P1prime_indices = which(all_orthologs$both_species_P4_P1prime == 1) 

cutsites$both_species_exact = NA
for (i in 1:nrow(cutsites)) {
  for (j in both_species_P4_P1prime_indices) {
    if(grepl(as.vector(all_orthologs$human_uniprot[j]), as.vector(cutsites$Entry.name[i]))) {
      cutsites$both_species_exact[i] = 1
    }   
  }
  print(i)
}  

# calculate number of unique proteins with cutsites$both_species_exact == 1
both_species_exact = unique(all_orthologs$human_uniprot[both_species_P4_P1prime_indices])
both_species_exact_total = length(both_species_exact)

## for conserved_shell 2: need to work out sensitivity threshold first - see later in script

# start sequence ranking
cutsites$fimo_rank = rank(-cutsites$score, ties.method = "min")

# extract known substrates from cutsites
fimo_ranking_malt1_substrates = cutsites[cutsites$Entry.name %in% known_substrates_entry_name,]
fimo_ranking_malt1_substrates_cutsites = 
  fimo_ranking_malt1_substrates[fimo_ranking_malt1_substrates$Matched.sequence %in% cutsite_coop$P4_P1prime,]
fimo_ranking_malt1_substrates_cutsites$sequence_name = as.character(fimo_ranking_malt1_substrates_cutsites$Entry.name)
fimo_ranking_malt1_substrates_cutsites$rank = fimo_ranking_malt1_substrates_cutsites$fimo_rank

#remove _HUMAN
fimo_ranking_malt1_substrates_cutsites$sequence_name = 
  sub("_HUMAN", "", fimo_ranking_malt1_substrates_cutsites$sequence_name)

# calculate sensitivity
number_of_cutsites = length(fimo_ranking_malt1_substrates_cutsites$Matched.sequence)

fimo_ranking_malt1_substrates_cutsites = fimo_ranking_malt1_substrates_cutsites[order(fimo_ranking_malt1_substrates_cutsites$rank),]

for(i in 1:number_of_cutsites) {
  fimo_ranking_malt1_substrates_cutsites$sensitivity[i] = i/number_of_cutsites
}

fimo_ranking_malt1_substrates_cutsites$sensitivity_min[1] = 0
fimo_ranking_malt1_substrates_cutsites$rank_min[1] = 0
for(i in 2:number_of_cutsites) {
  fimo_ranking_malt1_substrates_cutsites$sensitivity_min[i] = fimo_ranking_malt1_substrates_cutsites$sensitivity[i-1]
  fimo_ranking_malt1_substrates_cutsites$rank_min = fimo_ranking_malt1_substrates_cutsites$rank[i-1]
}

# use loess to determine FIMO ranks for sensitivity 0.5 and 0.8 of known MALT1 substrates
loess_fit <- loess(fimo_ranking_malt1_substrates_cutsites$rank ~ fimo_ranking_malt1_substrates_cutsites$sensitivity, fimo_ranking_malt1_substrates_cutsites)
sensitivity_output_0.5 = predict(loess_fit, newdata=0.5)
sensitivity_output_0.5 = floor(sensitivity_output_0.5)
sensitivity_output_0.8 = predict(loess_fit, newdata=0.8)
sensitivity_output_0.8 = floor(sensitivity_output_0.8)

#first make sure ranked by FIMO score, in case ordering changed
cutsites = cutsites[order(-cutsites$score),]

# get list of unique candidate cutsites below 0.8 sensitivity threshold
candidate_P4_P1prime = unique(cutsites$Matched.sequence[1:sensitivity_output_0.8])
# work out how many exact P4-P1prime matches in mouse proteome for calculation of score
all_orthologs$P4_P1prime_any_in_mouse = NA
all_orthologs$P4_P1prime_any_in_human = NA
for(i in 1:nrow(all_orthologs)) {
  for(j in 1:length(candidate_P4_P1prime)) {
    if(grepl(as.vector(candidate_P4_P1prime[j]), as.vector(all_orthologs$mouse_sequence[i]))) {
      all_orthologs$P4_P1prime_any_in_mouse[i] = 1
    }
    if(grepl(as.vector(candidate_P4_P1prime[j]), as.vector(all_orthologs$human_sequence[i]))) {
      all_orthologs$P4_P1prime_any_in_human[i] = 1
    }
  }
  print(i)
}

all_orthologs$both_species_P4_P1prime_any = NA
for(i in 1:nrow(all_orthologs)) {
  if(!is.na(all_orthologs$P4_P1prime_any_in_human[i])) {
    if(all_orthologs$P4_P1prime_any_in_human[i] == 1) {
      if(!is.na(all_orthologs$P4_P1prime_any_in_mouse[i])) {
        if(all_orthologs$P4_P1prime_any_in_mouse[i] == 1) {
          all_orthologs$both_species_P4_P1prime_any[i] = 1 } 
      } 
    }
  }
  print(i)
}

# append info onto cutsites
both_species_P4_P1prime_any_indices = which(all_orthologs$both_species_P4_P1prime_any == 1) 

cutsites$both_species_any = NA
for (i in 1:nrow(cutsites)) {
  for (j in both_species_P4_P1prime_any_indices) {
    if(grepl(as.vector(all_orthologs$human_uniprot[j]), as.vector(cutsites$Entry.name[i]))) {
      cutsites$both_species_any[i] = 1
    }   
  }
  print(i)
}  

# calculate number of unique proteins with cutsites$both_species_exact == 1
both_species_any = unique(all_orthologs$human_uniprot[both_species_P4_P1prime_any_indices])
both_species_any_total = length(both_species_any)

# Add conserved_scores onto cutsites
# Since at PROTEIN level, ensure for proteins with multiple candidate 
# cutsites, that conserved score only applied to that specific cutsite
# otherwise, bug with RC3H2
conserved_shell1_indices = which(cutsites$both_species_exact == TRUE & cutsites$Matched.sequence %in% cutsite_coop$P4_P1prime) 
conserved_shell2_indices = which(cutsites$both_species_any == TRUE)

# add in columns to output_uniprotKB, append interacts_shell info 
cutsites$conserved_shell1 = NA
cutsites$conserved_shell2 = NA
cutsites$conserved_shell1[conserved_shell1_indices] = 1
cutsites$conserved_shell2[conserved_shell2_indices] = 1

cutsites$conserved_shell1_score = NA
for(i in 1:nrow(cutsites)) {
  if (!is.na(cutsites$conserved_shell1[i])) {
    cutsites$conserved_shell1_score[i] = 1 / both_species_exact_total
  }
}

cutsites$conserved_shell2_score = NA
for(i in 1:nrow(cutsites)) {
  if (!is.na(cutsites$conserved_shell2[i])) {
    cutsites$conserved_shell2_score[i] = 1 / both_species_any_total
  }
}

# rank by FIMO score
cutsites$fimo_shell1 = NA
cutsites$fimo_shell2 = NA
cutsites$fimo_shell1_score = NA
cutsites$fimo_shell2_score = NA

cutsites$fimo_shell1[which(cutsites$fimo_rank < sensitivity_output_0.5)] = 1
cutsites$fimo_shell2[which(cutsites$fimo_rank < sensitivity_output_0.8)] = 1

unique_genes_fimo_shell1 = length(unique(cutsites$Entry.name[1:(sensitivity_output_0.5-1)]))
unique_genes_fimo_shell2 = length(unique(cutsites$Entry.name[1:(sensitivity_output_0.8-1)]))

cutsites$fimo_shell1_score = (cutsites$fimo_shell1 / unique_genes_fimo_shell1) * 
  (floor(number_known_substrates*0.5)/number_known_substrates)
cutsites$fimo_shell2_score = (cutsites$fimo_shell2 / unique_genes_fimo_shell2) * 
  (floor(number_known_substrates*0.8)/number_known_substrates)

### Sum sequence parameters
for(i in 1:nrow(cutsites)) {
  cutsites$sequence_combined_score[i] = sum(cutsites$exact_match_shell1_score[i], cutsites$exact_match_shell2_score[i], 
                                            cutsites$fimo_shell1_score[i], cutsites$fimo_shell2_score[i],
                                            cutsites$conserved_shell1_score[i],
                                            cutsites$conserved_shell2_score[i], na.rm=TRUE)
}

cutsites_sequence = cutsites

# order, so that 1st instance of each protein in cutsites_sequences is best sequence_rank 
cutsites_sequence = cutsites_sequence[order(cutsites_sequence$sequence_combined_score, decreasing = TRUE),]

# keep only best scoring sequence
cutsites_sequence = cutsites_sequence[match(unique(cutsites_sequence$Entry.name), cutsites_sequence$Entry.name),]

# rank by sequence_combined_score
cutsites_sequence$sequence_rank = rank(-cutsites_sequence$sequence_combined_score, ties.method = "min")

# convert cutsites to gene/protein level, rather than cleavage site level
cutsites = cutsites_sequence

# assign fimo rank (at gene level)
cutsites$fimo_rank = rank(-cutsites$score, ties.method = "min")

# load data.table package
library(data.table)
ppp = frank(cutsites, sequence_rank, fimo_rank, ties.method = 'min')
cutsites$sequence_rank_frank = ppp
cutsites$sequence_rank = cutsites$sequence_rank_frank

# apply unity-based normalisation to bio, sequence and fimo ranks
cutsites$bio_rank_normalised = (cutsites$bio_rank - min(cutsites$bio_rank)) / (max(cutsites$bio_rank) - min(cutsites$bio_rank))
cutsites$sequence_rank_normalised = (cutsites$sequence_rank - min(cutsites$sequence_rank)) / (max(cutsites$sequence_rank) - min(cutsites$sequence_rank))
cutsites$fimo_rank_normalised = (cutsites$fimo_rank - min(cutsites$fimo_rank)) / (max(cutsites$fimo_rank) - min(cutsites$fimo_rank))

# combine bio and sequence scores
cutsites$final_rank_normalised = cutsites$bio_rank_normalised * cutsites$sequence_rank_normalised

### Calculate GO-2-Substrates rank ####
cutsites = cutsites[order(cutsites$final_rank_normalised),]
cutsites_subset_unique_GO2Substrates_rank = cutsites[match(unique(cutsites$Gene_name), cutsites$Gene_name),]
cutsites_subset_unique_GO2Substrates_rank$gene_rank = rank(cutsites_subset_unique_GO2Substrates_rank$final_rank_normalised, ties.method = "min")
cutsites_subset_unique_GO2Substrates_rank$go2substrates_rank = cutsites_subset_unique_GO2Substrates_rank$gene_rank

write.csv(cutsites_subset_unique_GO2Substrates_rank, "Supplementary_Table_winnowing_mode.csv")


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

#### Calculate unique sequences and proteins for Fig 3c and 3d (Function Module 'bio'
#### and Sequence Module....note - FIMO must be calculated in Script A1)

# based on sequences rank (at cutsite level)
cutsites = cutsites[order(cutsites$sequence_rank),]
cutsites_subset_P4P1prime_sequence_rank = cutsites
cutsites_subset_P4P1prime_sequence_rank$gene_rank = rank(cutsites$sequence_rank, ties.method = "min")
max_sequence_rank_P4P1prime = max(which(cutsites_subset_P4P1prime_sequence_rank$literature_cleaved==1 &
                                          !is.na(cutsites_subset_P4P1prime_sequence_rank$exact_match_shell2_score)))
# obtain sequence rank at this criteria
temp_sequence = cutsites_subset_P4P1prime_sequence_rank$sequence_rank[max_sequence_rank_P4P1prime]
# subset those with sequence rank equal or better
max_sequence_rank_P4P1prime = max(which(cutsites_subset_P4P1prime_sequence_rank$sequence_rank == temp_sequence))
sequence_P4P1prime = cutsites_subset_P4P1prime_sequence_rank$Matched.sequence[1:max_sequence_rank_P4P1prime]
sequence_P4P1prime_unique = unique(sequence_P4P1prime)
# remove known cleavage sites to yield only CANDIDATE cutsites
sequence_P4P1prime_unique = sequence_P4P1prime_unique[! sequence_P4P1prime_unique %in% cutsite_coop$P4_P1prime]
write.csv(sequence_P4P1prime_unique, "unique_P4P1prime_sequence.csv")

# based on bio rank (at cutsite level)
cutsites = cutsites[order(cutsites$bio_rank),]
cutsites_subset_P4P1prime_bio_rank = cutsites
cutsites_subset_P4P1prime_bio_rank$gene_rank = rank(cutsites$bio_rank, ties.method = "min")
max_bio_rank_P4P1prime = max(which(cutsites_subset_P4P1prime_bio_rank$literature_cleaved==1 &
                                     !is.na(cutsites_subset_P4P1prime_bio_rank$exact_match_shell2_score)))
# obtain bio rank at this criteria
temp_bio = cutsites_subset_P4P1prime_bio_rank$bio_rank[max_bio_rank_P4P1prime]
# subset those with bio rank equal or better
max_bio_rank_P4P1prime = max(which(cutsites_subset_P4P1prime_bio_rank$bio_rank == temp_bio))
bio_P4P1prime = cutsites_subset_P4P1prime_bio_rank$Matched.sequence[1:max_bio_rank_P4P1prime]
bio_P4P1prime_unique = unique(bio_P4P1prime)
# remove known cleavage sites to yield only CANDIDATE cutsites
bio_P4P1prime_unique = bio_P4P1prime_unique[! bio_P4P1prime_unique %in% cutsite_coop$P4_P1prime]
write.csv(bio_P4P1prime_unique, "unique_P4P1prime_function.csv")



# based on sequences rank (at GENE / PROTEIN level)
cutsites = cutsites[order(cutsites$sequence_rank),]
cutsites_subset_P4P1prime_sequence_rank = cutsites
cutsites_subset_P4P1prime_sequence_rank$gene_rank = rank(cutsites$sequence_rank, ties.method = "min")
max_sequence_rank_P4P1prime = max(which(cutsites_subset_P4P1prime_sequence_rank$literature_cleaved==1 &
                                          !is.na(cutsites_subset_P4P1prime_sequence_rank$exact_match_shell2_score)))
# obtain sequence rank at this criteria
temp_sequence = cutsites_subset_P4P1prime_sequence_rank$sequence_rank[max_sequence_rank_P4P1prime]
# subset those with sequence rank equal or better
max_sequence_rank_P4P1prime = max(which(cutsites_subset_P4P1prime_sequence_rank$sequence_rank == temp_sequence))
sequence_P4P1prime = cutsites_subset_P4P1prime_sequence_rank$Entry.name[1:max_sequence_rank_P4P1prime]
sequence_P4P1prime_unique = unique(sequence_P4P1prime)
# remove known cleavage sites to yield only CANDIDATE cutsites
sequence_P4P1prime_unique = sequence_P4P1prime_unique[! sequence_P4P1prime_unique %in% known_substrates_entry_name]
write.csv(sequence_P4P1prime_unique, "unique_proteins_sequence.csv")


# based on bio rank (at GENE / PROTEIN level)
cutsites = cutsites[order(cutsites$bio_rank),]
cutsites_subset_P4P1prime_bio_rank = cutsites
cutsites_subset_P4P1prime_bio_rank$gene_rank = rank(cutsites$bio_rank, ties.method = "min")
max_bio_rank_P4P1prime = max(which(cutsites_subset_P4P1prime_bio_rank$literature_cleaved==1 &
                                     !is.na(cutsites_subset_P4P1prime_bio_rank$exact_match_shell2_score)))
# obtain bio rank at this criteria
temp_bio = cutsites_subset_P4P1prime_bio_rank$bio_rank[max_bio_rank_P4P1prime]
# subset those with bio rank equal or better
max_bio_rank_P4P1prime = max(which(cutsites_subset_P4P1prime_bio_rank$bio_rank == temp_bio))
bio_P4P1prime = cutsites_subset_P4P1prime_bio_rank$Entry.name[1:max_bio_rank_P4P1prime]
bio_P4P1prime_unique = unique(bio_P4P1prime)
# remove known cleavage sites to yield only CANDIDATE cutsites
bio_P4P1prime_unique = bio_P4P1prime_unique[! bio_P4P1prime_unique %in% known_substrates_entry_name]
write.csv(bio_P4P1prime_unique, "unique_proteins_function.csv")




