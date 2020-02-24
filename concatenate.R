## script: concatenate.R
## Concatenate all tables saved in QMAPP folder
library(tidyr)
library(dplyr)

# get file names 
files = list.files(path="QMAPP")
setwd("QMAPP")

exclude <- c("ABT199","Dexamethasone","Fingolimod","BMS345541","Ouabain", "Belinostat") 
fi = files[!grepl(paste(exclude, collapse="|"), files)]

# read all files into one dataframe
df = lapply(f, read.table, header=TRUE, sep=" ")
for (i in 1:length(df)){
  df[[i]] <- cbind(df[[i]], f[i])
  df[[i]]$dose_rep <- rownames(df[[i]])}

df <- do.call("rbind", df)
names(df)[14] <- "file"

# grab cell_line and drug columns
df <- tidyr::separate(data=df, col=file, into=c("cell_line","drug"), sep="\\.", extra="drop")

# grab dose and replicate
df <- tidyr::separate(data=df, col=dose_rep, into=c("dose","rep"))

# save (get rid of row names)
rownames(df) <- NULL
write.csv(df, file="QMAPP.concatenated.csv", row.names=FALSE)
