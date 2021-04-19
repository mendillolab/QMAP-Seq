# script: tidy.R
# tidy concatenated QMAP file into 3 different shapes/arrangements
# QMAP.concatenated.long // QMAP.doseResponse // QMAP.concatenated.AUC
# input: QMAP.concatenated.csv (built from concatenate.R)
# output: outpath (where to store results)

require(readr)
require(tidyr)
require(dplyr)

#### handle file paths #####
args <- commandArgs()
input <- sub('--input=', '', args[grep('--input=', args)])
output <- sub('--output=', '', args[grep('--output=', args)])

outLong <- paste0(output,"/QMAP.concatenated.long.csv")
outDR <- paste0(output,"/QMAP.concatenated.doseResponse.csv")
outAUC <- paste0(output,"/QMAP.concatenated.AUC.csv")

#### qmap long #####
qmap <- read_csv(input)
qmap.long <- qmap %>% tidyr::gather(sgRNA, survival, sgNT0.4:sgXBP1)
write.csv(qmap.long, file=outLong, row.names=FALSE)

#### qmap dose response ##### 
# load data and combine cell_line /drug columns
qmap <- qmap %>% unite(cell_line_drug, cell_line, drug, sep="_")

# convenience function: splitting guide columns into replicates
splitReps <- function(qmap, guide) {
  rep1 <-  qmap %>% select(guide, cell_line_drug, dose, rep) %>% filter(rep=="Rep1") %>% rename(rep1=guide)
  rep2 <-  qmap %>% select(guide, cell_line_drug, dose, rep) %>% filter(rep=="Rep2") %>% rename(rep2=guide)
  full <- full_join(rep1, rep2, by=c("cell_line_drug", "dose")) %>% select(-rep.x, -rep.y)
  colnames(full) <- c(paste0(guide,"_rep1"), "cell_line_drug", "dose", paste0(guide, "_rep2"))
  return(full)
}

# for each guide, combine into aggregate df
guides <- c("sgATF3","sgATF4","sgATF6","sgATG7","sgERN1","sgHSF1","sgHSF2","sgKEAP1","sgNFE2L2","sgSLC35F2","sgXBP1")
qmap.DR <- splitReps(qmap, "sgNT0.4") # initialize w/ non-target data
for (guide in guides) {
  split.df <- splitReps(qmap, guide)
  qmap.DR <- left_join(qmap.DR, split.df, by=c("cell_line_drug", "dose"))
}

# re-order columns & save
column_order <- c("cell_line_drug","dose", "sgNT0.4_rep1", "sgNT0.4_rep2", 
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
qmap.DR <- qmap.DR[,column_order]
write_csv(qmap.DR, path=outDR)

#### qmap AUC ##### 
# AUC = 0.5*Dose1 + Dose2 + Dose3 + 0.5*Dose4
qmap.AUC <- qmap.long 
qmap.AUC.1 <- qmap.AUC %>% filter(rep=="Rep1") %>% spread(dose, survival) %>% mutate(AUC.1=0.5*Dose1+Dose2+Dose3+0.5*Dose4)
qmap.AUC.2 <- qmap.AUC %>% filter(rep=="Rep2") %>% spread(dose, survival) %>% mutate(AUC.2=0.5*Dose1+Dose2+Dose3+0.5*Dose4)
qmap.AUC <- cbind(qmap.AUC.1, qmap.AUC.2[,"AUC.2"]) %>% select(-rep, -(Dose1:Dose4))
write_csv(qmap.AUC, path=outAUC)
print("Finished tidy.R")



