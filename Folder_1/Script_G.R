##### Script G #####
##### R script to generate MALT1 PSSM and reproduce Fig 1e in GO-2-Substrates paper #####

# load P4-P1' sequence info for each knwon MALT1 cleavage site
MALT1_Cleavages = read.csv("cleavage_sites_for_PSSM.csv", header=TRUE)

# first, derive P4...P3...etc for each P4P1prime entry
test_temporary_frame = MALT1_Cleavages

Ppositions = c("P4", "P3", "P2", "P1", 
               "P1prime")
Ppositions = as.vector(Ppositions)
amino_acids = c("A", "C", "D", "E", "F","G","H","I","K","L",
                "M","N","P","Q","R","S","T","V","W","Y")
amino_acids = as.character(amino_acids)
temp_frame = data.frame(matrix(, nrow=5, ncol=20))   
colnames(temp_frame) = amino_acids
rownames(temp_frame) = Ppositions

# then iterate over P positions, calculate frequency for each position of cleavage site (P4 to P1') for each amino acid
P_tmp_prop_vector_combined = NA
for (j in (Ppositions)) {
  P_tmp = table(test_temporary_frame[,j])
  P_tmp_prop = P_tmp/sum(P_tmp)
  # Create vector for each position, assigning frequencies to each of 20 natural amino acids, in alphabetical order  
  P_tmp_prop_vector = c(P_tmp_prop["A"], P_tmp_prop["C"], P_tmp_prop["D"], P_tmp_prop["E"], P_tmp_prop["F"],
                        P_tmp_prop["G"], P_tmp_prop["H"], P_tmp_prop["I"], P_tmp_prop["K"], P_tmp_prop["L"],
                        P_tmp_prop["M"], P_tmp_prop["N"], P_tmp_prop["P"], P_tmp_prop["Q"], P_tmp_prop["R"],
                        P_tmp_prop["S"], P_tmp_prop["T"], P_tmp_prop["V"], P_tmp_prop["W"], P_tmp_prop["Y"])
  for(z in 1:length(P_tmp_prop_vector)) {
    if(is.na(P_tmp_prop_vector[z])) {
      P_tmp_prop_vector[z] = 0
    }
  }
  temp_frame[j,] = P_tmp_prop_vector
  temp_frame_combined = temp_frame
  rm(P_tmp)
  #P_tmp_prop_vector_combined =  subset(P_tmp_prop_vector_combined, select = -r[1])
}
temp_protease = "total"
temp_index = "total"
total_PSSM = temp_frame_combined

write.csv(total_PSSM, "Fig_1e.csv")
