####################################################
### QMAPP reformatting for PRISM plots - 10/8/19 ###

# load data and combine cell_line /drug columns
qmapp <- read_csv("QMAPP/QMAPP.concatenated.csv")
qmapp <- qmapp %>% unite(cell_line_drug, cell_line, drug, sep="_")
qmapp <- qmapp %>% mutate(logDrugConcentration=log10(Concentration)) %>% 
                   select(-Concentration)

# convenience function: splitting guide columns into replicates
splitReps <- function(qmapp, guide) {
  rep1 <-  qmapp %>% select(logDrugConcentration, guide, cell_line_drug, dose, rep) %>% filter(rep=="Rep1") %>% rename(rep1=guide)
  rep2 <-  qmapp %>% select(logDrugConcentration, guide, cell_line_drug, dose, rep) %>% filter(rep=="Rep2") %>% rename(rep2=guide)
  full <- full_join(rep1, rep2, by=c("cell_line_drug", "dose", "logDrugConcentration")) %>% select(-rep.x, -rep.y)
  colnames(full) <- c("logDrugConcentration", paste0(guide,"_rep1"), "cell_line_drug", "dose", paste0(guide, "_rep2"))
  return(full)
}

# for each guide, combine into aggregate df
guides <- c("sgATF3","sgATF4","sgATF6","sgATG7","sgERN1","sgHSF1","sgHSF2","sgKEAP1","sgNFE2L2","sgSLC35F2","sgXBP1")
qmapp2 <- splitReps(qmapp, "sgNT0.4") # initialize w/ non-target data
for (guide in guides) {
  split.df <- splitReps(qmapp, guide)
  qmapp2 <- left_join(qmapp2, split.df, by=c("cell_line_drug", "dose", "logDrugConcentration"))
}

# re-order columns & save
column_order <- c("cell_line_drug","logDrugConcentration","dose", "sgNT0.4_rep1", "sgNT0.4_rep2", 
                                                           "sgATF3_rep1", "sgATF3_rep2",
                                                           "sgATF4_rep1", "sgATF4_rep2",
                                                           "sgATF6_rep1", "sgATF6_rep2",
                                                           "sgATG7_rep1", "sgATG7_rep2",
                                                           "sgERN1_rep1", "sgERN1_rep2",
                                                           "sgHSF1_rep1", "sgHSF1_rep2",
                                                           "sgHSF2_rep1", "sgHSF2_rep2",
                                                           "sgKEAP1_rep1", "sgKEAP1_rep2",
                                                           "sgNFE2L2_rep1", "sgNFE2L2_rep2",
                                                           "sgSLC35F2_rep1", "sgSLC35F2_rep2",
                                                           "sgXBP1_rep1", "sgXBP1_rep2") 
qmapp2 <- qmapp2[,column_order]
write_csv(qmapp2, path="QMAPP.concatenated.forDoseResponse.csv")
