
### Work out orthologues and append to filter_pool_A_B1_B2_C1_C2_all
# import orthologues
all_orthologs = read.csv(file="all_orthologs_161120.csv", header=TRUE)

# import uniprotKB results for all human ensemble genes matching to swissprot protein
uniprot_reviewed_ensemble_human = read.csv(file="uniprot_reviewed_human_ensemble_161120.csv", header=TRUE)
# import uniprotKB results for all human ensemble genes matching to unreviewed uniprot protein
uniprot_unreviewed_ensemble_human = read.csv(file="uniprot_unreviewed_human_ensemble_161120.csv", header=TRUE)
# remove all fragments from unreviewed entries, to avoid matching on incorrect entry
uniprot_unreviewed_ensemble_human = uniprot_unreviewed_ensemble_human %>% filter(str_detect(uniprot_unreviewed_ensemble_human$Protein.names, "(Fragment)", negate = TRUE))
# import uniprotKB results for all mouse ensemble genes matching to swissprot protein
uniprot_reviewed_ensemble_mouse = read.csv(file="uniprot_reviewed_mouse_ensemble_161120.csv", header=TRUE)
# import uniprotKB results for all mouse ensemble genes matching to unreviewed uniprot protein
uniprot_unreviewed_ensemble_mouse = read.csv(file="uniprot_unreviewed_mouse_ensemble_161120.csv", header=TRUE)
# remove all fragments from unreviewed entries, to avoid matching on incorrect entry
uniprot_unreviewed_ensemble_mouse = uniprot_unreviewed_ensemble_mouse %>% filter(str_detect(uniprot_unreviewed_ensemble_mouse$Protein.names, "(Fragment)", negate = TRUE))


# 1st, match ensemble ids to Uniprot IDs of reviewed (Swissprot) proteins
reviewed_human_indices = which(all_orthologs$human %in% uniprot_reviewed_ensemble_human$Ensemble)
reviewed_human_uniprot_indices = which(uniprot_reviewed_ensemble_human$Ensemble %in% all_orthologs$human)

all_orthologs$human_uniprot = NA
all_orthologs$human_sequence = NA
temp_vector = NA
temp_sequence = NA
for (i in reviewed_human_indices) {
  for (j in reviewed_human_uniprot_indices) {
    if(identical(as.vector(all_orthologs$human[i]), as.vector(uniprot_reviewed_ensemble_human$Ensemble[j]))) {
      temp_vector = as.character(uniprot_reviewed_ensemble_human$Entry.name[j])
      temp_sequence = as.character(uniprot_reviewed_ensemble_human$Sequence[j])
      if(!is.na(temp_vector)) break }
  }
  all_orthologs$human_uniprot[i] = temp_vector
  all_orthologs$human_sequence[i] = temp_sequence
  temp_vector = NA
  temp_sequence = NA
  print(i)
}

# next, if all_orthologs$human_uniprot still empty, match against unreviewed uniprot IDs
unreviewed_human_indices = which(all_orthologs$human %in% uniprot_unreviewed_ensemble_human$Ensemble)
unreviewed_human_uniprot_indices = which(uniprot_unreviewed_ensemble_human$Ensemble %in% all_orthologs$human)

temp_vector = NA
for (i in unreviewed_human_indices) {
  if(is.na(all_orthologs$human_uniprot[i])) {
    for (j in unreviewed_human_uniprot_indices) {
      if(identical(as.vector(all_orthologs$human[i]), as.vector(uniprot_unreviewed_ensemble_human$Ensemble[j]))) {
        temp_vector = as.character(uniprot_unreviewed_ensemble_human$Entry.name[j]) 
        temp_sequence = as.character(uniprot_unreviewed_ensemble_human$Sequence[j])
        if(!is.na(temp_vector)) break }
    }}
  if(!is.na(temp_vector)) {
    all_orthologs$human_uniprot[i] = temp_vector
    all_orthologs$human_sequence[i] = temp_sequence
  }
  temp_vector = NA
  temp_sequence = NA
  print(i)
}


all_orthologs_backup3 = all_orthologs


# import uniprotKB results for all mouse ensemble genes matching to swissprot protein
# get row numbers of all_orthologs present in uniprot_reviewed_ensemble_mouse
reviewed_mouse_indices = which(all_orthologs$mouse %in% uniprot_reviewed_ensemble_mouse$Ensemble)
reviewed_mouse_uniprot_indices = which(uniprot_reviewed_ensemble_mouse$Ensemble %in% all_orthologs$mouse)
# 1st, match ensemble ids to Uniprot IDs of reviewed (Swissprot) proteins
all_orthologs$mouse_uniprot = NA
all_orthologs$mouse_sequence = NA
temp_vector = NA
temp_sequence = NA
for (i in reviewed_mouse_indices) {
  for (j in reviewed_mouse_uniprot_indices) {
    if(identical(as.vector(all_orthologs$mouse[i]), as.vector(uniprot_reviewed_ensemble_mouse$Ensemble[j]))) {
      temp_vector = as.character(uniprot_reviewed_ensemble_mouse$Entry.name[j])
      temp_sequence = as.character(uniprot_reviewed_ensemble_mouse$Sequence[j])
      if(!is.na(temp_vector)) break }
  }
  all_orthologs$mouse_uniprot[i] = temp_vector
  all_orthologs$mouse_sequence[i] = temp_sequence
  temp_vector = NA
  temp_sequence = NA
  print(i)
}

all_orthologs_backup4 = all_orthologs

# get row numbers of all_orthologs present in uniprot_reviewed_ensemble_mouse
unreviewed_mouse_indices = which(all_orthologs$mouse %in% uniprot_unreviewed_ensemble_mouse$Ensemble)
unreviewed_mouse_uniprot_indices = which(uniprot_unreviewed_ensemble_mouse$Ensemble %in% all_orthologs$mouse)
# next, if all_orthologs$mouse_uniprot still empty, match against unreviewed uniprot IDs
temp_vector = NA
temp_sequence = NA
for (i in unreviewed_mouse_indices) {
  for (j in unreviewed_mouse_uniprot_indices) {
    if(is.na(all_orthologs$mouse_uniprot[i])) {
      if(identical(as.vector(all_orthologs$mouse[i]), as.vector(uniprot_unreviewed_ensemble_mouse$Ensemble[j]))) {
        temp_vector = as.character(uniprot_unreviewed_ensemble_mouse$Entry.name[j])
        temp_sequence = as.character(uniprot_unreviewed_ensemble_mouse$Sequence[j])
        if(!is.na(temp_vector)) break }
    }}
  if(!is.na(temp_vector)) { all_orthologs$mouse_uniprot[i] = temp_vector }
  if(!is.na(temp_sequence)) { all_orthologs$mouse_sequence[i] = temp_sequence }
  temp_vector = NA
  temp_sequence = NA
  print(i)
}

write.csv(all_orthologs, "orthologs.csv")
