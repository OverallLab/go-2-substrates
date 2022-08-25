##### Script H #####
##### R script to generate PPI input used for Cytoscape visualisation, using BioGrid data #####

# script to count number of instances of evidence for each PPI for each new MALT1 substrate
file_inputs = c("BIOGRID-GENE-107293-4.3.196_CASP10_human_only.tab3.txt",
                "BIOGRID-GENE-115328-4.3.196_TANK_human_only.tab3.txt",
                "BIOGRID-GENE-116527-4.3.196_CILK1.tab3.txt",
                "BIOGRID-GENE-129216-4.3.196_TAB3_human_only.tab3.txt",
                "BIOGRID-GENE-131004-4.3.196_ZC3H12D.tab3.txt",
                "BIOGRID-GENE-131074-4.3.196_ZC3H12B.tab3.txt",
                "BIOGRID-GENE-132364-4.3.196_ILDR2_human_only.tab3.txt")

for(j in 1:length(file_inputs)) {
  # import biogrid file
  biogrid_download = read.delim(file_inputs[j], header=TRUE)
  # work out sum and product of IDs, to be used as unique IDs in combination
  biogrid_download$BioGRID_ID_sumAB = as.numeric(biogrid_download$BioGRID.ID.Interactor.A) + as.numeric(biogrid_download$BioGRID.ID.Interactor.B)
  biogrid_download$BioGRID_ID_productAB = as.numeric(biogrid_download$BioGRID.ID.Interactor.A) * as.numeric(biogrid_download$BioGRID.ID.Interactor.B, na.rm=TRUE) 
  
  # filter all rows with the same sumAB and productAB (i.e. are the same PPI)
  # count total number of evidence for interaction, number of publications, unique experimental systems
  for(i in 1:nrow(biogrid_download)) {
    temp = biogrid_download %>% filter(BioGRID_ID_sumAB[i] == BioGRID_ID_sumAB &
                                         BioGRID_ID_productAB[i] == BioGRID_ID_productAB)
    temp_evidence_count = nrow(temp)
    temp_experimental_system_count = n_distinct(temp$Experimental.System)
    temp_publication_count = n_distinct(temp$Publication.Source)
    biogrid_download$Evidence_count[i] = temp_evidence_count
    biogrid_download$experimental_system_count[i] = temp_experimental_system_count
    biogrid_download$publication_count[i] = temp_publication_count
  }
  # subset only human PPIs, physical interactors
  biogrid_download_subset = biogrid_download %>% filter(Organism.ID.Interactor.A == "9606" &
                                                          Organism.ID.Interactor.B == "9606" )
  biogrid_download_subset = biogrid_download_subset %>% filter(Experimental.System.Type == "physical")
  biogrid_download_subset = biogrid_download_subset %>% select(Official.Symbol.Interactor.A, Official.Symbol.Interactor.B,
                                                               Experimental.System, Experimental.System.Type,
                                                               Author, Publication.Source, Throughput,
                                                               BioGRID_ID_sumAB, BioGRID_ID_productAB, 
                                                               Organism.ID.Interactor.A, Organism.ID.Interactor.B,
                                                               BioGRID.ID.Interactor.A, BioGRID.ID.Interactor.B,
                                                               Evidence_count, experimental_system_count,
                                                               publication_count)
  write.csv(biogrid_download_subset, paste0("150421_Evidence_count",file_inputs[j],".csv"), quote = FALSE)
  print(j)
}