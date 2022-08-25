##### Script D #####
##### R script to reproduce Figures 1b, 1c, 1f, Supplementary 8h, in GO-2-Substrates paper #####


#### Plot Fig 1b
human_pre_study = c("PLVPRGGGTP","QLIPRGTDPS","QMVPRGSQLY","QLISRSTDST",
                    "CHCSRTPDAF","KLCCRATGHP","SCLSRGAHEK","AFMSRGVGDK",
                    "RLVSRGAASL","PLRSRTVSRQ","LGASRGEAYE","TLQPRGPLEP",
                    "SPDSRASSLS","LFKSKGNYDE","NFVSRGASSH")
human_pre_study_logo = ggseqlogo(as.character(human_pre_study))
human_pre_study_logo +
  ggsave("Fig_1b.pdf")

#### Plot Fig 1c
mouse_pre_study = c("PLVPRGGSTP","QLIPRGTDPS","QMVPRGSQLY",
                    "QLIPRGTDSA","CHCSRTPHTF","KLCCRATGHP",
                    "SCPSRGALEK","AFMSRGVGDK","RLVPRGPASL",
                    "PLRSRALSRQ","TLQSRGPLEP","LFKSKGNYDE",
                    "SFVSRGASGH")

mouse_pre_study_logo = ggseqlogo(as.character(mouse_pre_study))
mouse_pre_study_logo +
  ggsave("Fig_1c.pdf")

#### Plot Supplementary Fig. 8h
library(ggseqlogo)
negative_from_screen = c("ALLSRGGVYA", "LLVPRAVGKI", "VLVSRGGSQS", "LLVSRGQRRL", "QFCSRGFREK",
                         "CLVPRGRRLG", "VLCSRGTRAK", "NMASRGNINV", "PEPSRGGVSV", "GPVSRGRGLL",
                         "RLQSRTEDSD", "LFKSKGNYDE", "LGASRGSSVE", "GFVPRGAFGK","FLVPRGTKME",
                         "LLQPRGGSQE", "SPRPRAFDGA", "CLRSRTPYHV", "APVSRGRDGY", "DPDSRGCCQE",
                         "ILDPKGSLLG", "SLRSRAKGRL", "RLVPRTGRGA", "DLKPRGLTGM", "HLQPRGQPAP",
                         "PPVPRGRGVG", "ELVPRSILAM", "IPLPRGQTEK", "SPQPRSTPRQ", "PPRSRGPPRG",
                         "LLKSKGTCQS", "VGPPRGALVR", "GPPSRGGHMD", "RGAPRGGGRG")
ggseqlogo(negative_from_screen)
ggsave(width = 16, height = 12, units = "cm", filename ="Supplementary_Fig8h.pdf")

#### Plot Fig. 1f

# rank scores from fimo
fimo_ranking = read.csv("fimo_sensitivity_1.csv", header = TRUE)
fimo_ranking$rank = NA
fimo_ranking$rank = rank(-fimo_ranking$score, ties.method = "min")

# define substrate IDs
known_substrates_entry_name = 
  c("ZC12A_HUMAN", "RC3H1_HUMAN", "HOIL1_HUMAN", "CYLD_HUMAN", "BCL10_HUMAN", 
    "RC3H2_HUMAN", "M3K14_HUMAN", "LIMA1_HUMAN", "MALT1_HUMAN", "TNAP3_HUMAN",
    "RELB_HUMAN", "N4BP1_HUMAN")
# extract known substrates from cutsites
cutsites[cutsites$Entry.name %in% known_substrates_entry_name,]
# load known cleavage sites
cutsite_coop = read.csv("cleavage_sites.csv",header=TRUE)
#subset rows in fimo_ranking with IDs that are equal to MALT1 substrates
fimo_ranking_malt1_substrates = fimo_ranking[fimo_ranking$sequence_name %in% known_substrates_entry_name,]
fimo_ranking_malt1_substrates_cutsites = 
  fimo_ranking_malt1_substrates[fimo_ranking_malt1_substrates$matched_sequence %in% cutsite_coop$P4_P1prime,]
fimo_ranking_malt1_substrates_cutsites$sequence_name = as.character(fimo_ranking_malt1_substrates_cutsites$sequence_name)
# remove _HUMAN from IDs
fimo_ranking_malt1_substrates_cutsites$sequence_name = 
  sub("_HUMAN", "", fimo_ranking_malt1_substrates_cutsites$sequence_name)

# Prepare data for plotting
fimo_ranking_malt1_substrates_cutsites$gene = c("RELB","N4BP1","ZC3H12A", "RC3H1-1", "CYLD", "RC3H1-2",
                                                "HOIL1",
                                                "MALT1-1", "MALT1-2", "MAP3K14", "BCL10",
                                                "RC3H2", "LIMA1-1", "TNFAIP3", "LIMA-2")

fimo_ranking_malt1_substrates_cutsites$gene_common = c("RELB","N4BP1","ZC3H12A", "RC3H1", "CYLD", "RC3H1",
                                                       "HOIL1",
                                                       "MALT1", "MALT1", "MAP3K14", "BCL10",
                                                       "RC3H2", "LIMA1", "TNFAIP3", "LIMA1")

fimo_ranking_malt1_substrates_cutsites$autocleavage = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                                        FALSE,
                                                        TRUE, TRUE, FALSE, FALSE,
                                                        FALSE, FALSE, FALSE, FALSE)

fimo_ranking_malt1_substrates_cutsites$cleaved_by = c("CBM","CBM","CBM", "CBM", "CBM", "CBM",
                                                      "CBM",
                                                      "CBM", "CBM", "cIAP2-MALT1", "CBM",
                                                      "CBM", "cIAP2-MALT1", "CBM", "cIAP2-MALT1")

fimo_ranking_malt1_substrates_cutsites$cutsite_number = c(1,1,1,1,1,2,1,1,2,1,1,1,1,1,2)

number_of_cutsites = length(fimo_ranking_malt1_substrates_cutsites$matched_sequence)

# Work out sensitivity and fimo rank at each cleavage site
for(i in 1:number_of_cutsites) {
  fimo_ranking_malt1_substrates_cutsites$sensitivity[i] = i/number_of_cutsites
}

fimo_ranking_malt1_substrates_cutsites$sensitivity_min[1] = 0
fimo_ranking_malt1_substrates_cutsites$rank_min[1] = 0

for(i in 2:number_of_cutsites) {
  fimo_ranking_malt1_substrates_cutsites$sensitivity_min[i] = fimo_ranking_malt1_substrates_cutsites$sensitivity[i-1]
  fimo_ranking_malt1_substrates_cutsites$rank_min = fimo_ranking_malt1_substrates_cutsites$rank[i-1]
}

##### generate Fig. 1f
ggplot(fimo_ranking_malt1_substrates_cutsites, aes(y=rank, x=sensitivity)) +
  scale_y_reverse(breaks = c(1, 1000, 2000, 3000, 4000, 5000, 6000), 
                  labels = c(1, 1000, 2000, 3000, 4000, 5000, 6000),
                  expand = c(0.02,0)) +
  coord_cartesian(ylim = c(6000, 1)) +
  scale_x_continuous(position = 'top', labels = c("0","0.25","0.5","0.75","1"), 
                     breaks = seq(0, 1, by = 0.25), limits = c(0,1), expand = c(0,0.019)) +
  geom_point(aes(color = factor(cleaved_by)),size=3) +
  stat_smooth(geom = "line",alpha = 0.3, color = "black", span = 0.7) +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x=element_text(size = 6, face = "bold"))+
  theme(axis.title.y=element_text(size = 6, face = "bold"))+
  theme(axis.text.x=element_text(size = 6, face = "bold")) +
  theme(axis.text.y=element_text(size = 6, face = "bold")) +
  theme(axis.line = element_line(color = "black",size = 0.5)) +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(size = 6, face = "bold")) +
  theme(legend.title = element_blank()) +
  theme(panel.background = element_blank()) +
  labs(x="Sensitivity", y="FIMO cutsite ranking", col="cleaved_by")
ggsave("Fig_1f.pdf", width = 90, height = 85, units = "mm")

