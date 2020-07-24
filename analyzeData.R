library(gplots)
library(drc)

args <- commandArgs()

pattern <-sub('--pattern=', '', args[grep('--pattern=', args)])
directory <- sub('--directory=', '', args[grep('--directory=', args)])
doselist <- sub('--doselist=','',args[grep('--doselist=', args)])
spikein <- sub('--spikein=','',args[grep('--spikein=', args)])
intermediary <- sub('--intermediary=','',args[grep('--intermediary=', args)])
if (identical(intermediary,character(0))){
  intermediary<- 0
}
print(intermediary)

doses<-read.table(doselist,header=TRUE,sep="\t",row.names=1)
doses
spikein<-read.table(spikein,header=TRUE,sep="\t")
spikein$logCellNumber <- log2(spikein$CellNumber)
controlCell <- toString(spikein$CellLine[1])
controlGuides <- spikein$sgRNA
print(paste("Control cell is",controlCell))
print("Control guides are")
print(controlGuides)
rownames(spikein)<-spikein$sgRNA
print(spikein)
spikein
numdoses <- ncol(doses)
print(paste("Num doses:",numdoses))
numreps<- 2
replabels<-c("Rep1","Rep2")

setwd(directory)
fileList <- list.files(pattern = pattern)

# First read in all of the data into a large matrix
for (file in fileList){
  fileLabel <- gsub(".labeledCounts.txt","",file)
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
     dataset<- read.table(file, header=FALSE, sep="\t")
     colnames(dataset)<-c("CellLine_Barcode","sgRNA_Barcode","CellLine","sgRNA",fileLabel)
  }    

  # if the merged dataset does exist, append to it
  else if (exists("dataset")){
     temp_dataset <-read.table(file, header=FALSE, sep="\t")
     dataset<-cbind(dataset,temp_dataset[,5])
     colnames(dataset)[ncol(dataset)] <- fileLabel
     rm(temp_dataset)
  }
}

# Fix the drug names and remove sample numbers, if necessary.
colnames(dataset)<-gsub("Taxmoxifen_Citrate","TamoxifenCitrate",colnames(dataset))	
colnames(dataset)<-gsub("Tamoxifen_Citrate","TamoxifenCitrate",colnames(dataset))
colnames(dataset)<-gsub("YM_155","YM155",colnames(dataset))
colnames(dataset)<-gsub("_S[0-9]+$","",colnames(dataset))

# Setup rownames with cellLine and sgRNA
rownames(dataset)<-paste(dataset$CellLine,dataset$sgRNA,sep=".")

# Extract just the read numbers.
justNumbers<-dataset[,5:ncol(dataset)]
colnames(justNumbers)<-gsub("SB_ZSHMB_sgPool_","",colnames(justNumbers))

# Remove extra sample information after Drug, to pull out the list of drugs
drugs <- gsub('_Dose[0-9]+_Rep[0-9]+$', '', colnames(justNumbers))
drugs <- unique(drugs)
print("Drugs are:")
print(drugs)

# Pull out DMSO samples
dmsoSamples <- drugs[grep("DMSO",drugs)]
# Remove DMSO from the drugs list
drugs <- drugs[grep("DMSO",drugs,invert=TRUE)]
numdrugs<-length(drugs)
print("nonDMSO drugs:")
print(drugs)
print("DMSO samples:")
print(dmsoSamples)

cellLines <- unique(dataset$CellLine)
cellLines<-cellLines[grep(controlCell,cellLines,invert=TRUE)]
print("CellLines:")
print(cellLines)

for (cellLine in cellLines){
   for (drug in drugs){
      for (dose in colnames(doses)){
	 for (rep in replabels){
	    sgRNAandCell0 <- paste(cellLine,"sgNT0",sep=".")
	    sgRNAandCell1 <- paste(cellLine,"sgNT1",sep=".")
	    sgRNAandCell2 <- paste(cellLine,"sgNT2",sep=".")
	    sgRNAandCell3 <- paste(cellLine,"sgNT3",sep=".")
	    sgRNAandCell4 <- paste(cellLine,"sgNT4",sep=".")
	    newlabel <- paste(cellLine,"sgNT0-4",sep=".")
#	    print(newlabel)
            colLabel <- paste(drug,dose,rep,sep="_")
#	    print(colLabel)
	    justNumbers[newlabel,colLabel]<-(justNumbers[sgRNAandCell0,colLabel]+justNumbers[sgRNAandCell1,colLabel]+justNumbers[sgRNAandCell2,colLabel]+justNumbers[sgRNAandCell3,colLabel]+justNumbers[sgRNAandCell4,colLabel])
	 }
      }
   }
}

print("myCols")
myCols <- colnames(justNumbers)
print(myCols)
print("Adding DMSO sample control numbers")
# Now for the DMSO samples
for (cellLine in cellLines){
   for (dmso in dmsoSamples){
       sgRNAandCell0 <- paste(cellLine,"sgNT0",sep=".")
       sgRNAandCell1 <- paste(cellLine,"sgNT1",sep=".")
       sgRNAandCell2 <- paste(cellLine,"sgNT2",sep=".")
       sgRNAandCell3 <- paste(cellLine,"sgNT3",sep=".")
       sgRNAandCell4 <- paste(cellLine,"sgNT4",sep=".")
       newlabel <- paste(cellLine,"sgNT0-4",sep=".")
       colLabel <- colnames(justNumbers)[grep(dmso,colnames(justNumbers))]
       justNumbers[newlabel,colLabel]<-(justNumbers[sgRNAandCell0,colLabel]+justNumbers[sgRNAandCell1,colLabel]+justNumbers[sgRNAandCell2,colLabel]+justNumbers[sgRNAandCell3,colLabel]+justNumbers[sgRNAandCell4,colLabel])
   }
}

print("justNumbers")
print(head(justNumbers))

# Take the log
print("Take the log")
logNumbers<-log2(justNumbers)
print("logNumbers")
print(head(logNumbers))

# Print intermediary numbers, for debugging.
if (intermediary == 1){
    dir.create("intermediary/")
   for (cellLine in cellLines){
      for (drug in drugs){
         subsetJustNumbers<-justNumbers[grep(cellLine,rownames(justNumbers)),grep(drug,colnames(logNumbers))]
         write.table(subsetJustNumbers,paste("intermediary/",cellLine,".",drug,".readCounts.txt",sep=""))
       	 subsetLogNumbers<-logNumbers[grep(cellLine,rownames(logNumbers)),grep(drug,colnames(logNumbers))]
      	 write.table(subsetLogNumbers,paste("intermediary/",cellLine,".",drug,".log2readCounts.txt",sep=""))
      }
      subsetJustNumbers<-justNumbers[grep(cellLine,rownames(justNumbers)),grep("DMSO",colnames(logNumbers))]
      write.table(subsetJustNumbers,paste("intermediary/",cellLine,".DMSOsamples.readCounts.txt",sep=""))
      subsetLogNumbers<-logNumbers[grep(cellLine,rownames(logNumbers)),grep("DMSO",colnames(logNumbers))]
      write.table(subsetLogNumbers,paste("intermediary/",cellLine,".DMSOsamples.log2readCounts.txt",sep=""))
   }
}

## DMSO samples still need to be worked through cell number interpolation

print("Check Spike-in")
# Use the spike in information to figure out how many cells were spiked in
print("Control Cell")
print(controlCell)

print(cellLines)
print(rownames(logNumbers))
spikeInData <-logNumbers[grepl(controlCell,rownames(logNumbers)),]
print(spikeInData)
rownames(spikeInData) <- gsub(controlCell,'',rownames(spikeInData))
rownames(spikeInData) <- gsub('^\\.','',rownames(spikeInData))
print(spikeInData)
numdrugs <- length(drugs)
numSamples <- numreps*numdrugs*numdoses
print(numdrugs)
print(numSamples)
interpolationRows <- colnames(logNumbers)
interpolationRows <- interpolationRows[grep("DMSO",interpolationRows,invert=TRUE)]
interpolationInfo<-data.frame(matrix(ncol=2,nrow=numSamples))
colnames(interpolationInfo)<-c("slope","intercept")
rownames(interpolationInfo)<-interpolationRows

print("drugs")
print(drugs)
# Use the spikein data to infer how to interpolate from read numbers to cell
# numbers for each well in the experiment.

problemDrugs <- c("Dexamethasone","ABT199","Fingolimod","BMS345541","Ouabain","Belinostat")
for (drug in problemDrugs){
    drugs <- drugs[grep(drug,drugs,invert=TRUE)]
}
#testDrugs <-c("FourHydroxytamoxifen","Lapatinib","YM155","Ganetespib")
#testDrugs <-c("YM155","Carfilzomib")
#drugs <-testDrugs

dir.create("plots")
for (drug in drugs){
   print(paste("Drug:", drug))
   spikeInDrug <- spikeInData[,grepl(drug,colnames(spikeInData))]
   spikeInDrug <- spikeInDrug[order(rownames(spikeInDrug)),]
   spikeInDrugControls <- spikeInDrug[grepl("^sgNT",rownames(spikeInDrug)),]
   spikeInDrugControls <- spikeInDrugControls[6:10,]
   if(identical(drug,drugs[1])){
      print("Rownames of selected Controls:")
      print(toString(rownames(spikeInDrugControls)))
      print("List of Controls from spike in file:")
      print(toString(controlGuides))
      print("If the two lists above don't match, this will lead to incorrect assumptions.")
   }
   pdf(paste("plots/",controlCell,".sgNT.",drug,".plot.pdf",sep=""))
   par(mar=c(12.1,4.1,2.1,6))
       matplot(spikein$logCellNumber,spikeInDrugControls,ylab="log2(readCounts)",xlab="log2(cellNumber)",pch=1,main=paste(controlCell,"sgNT",drug))
   nc<-ncol(spikeInDrugControls)
   print(spikeInDrugControls)
   for (i in 1:nc){
      sample <- colnames(spikeInDrugControls)[i]
      print(sample)
      lmSpikeIn <- lm(spikeInDrugControls[,i] ~ spikein$logCellNumber)
      interpolationInfo[sample,"slope"] <-coef(lmSpikeIn)["spikein$logCellNumber"]
      interpolationInfo[sample,"intercept"] <- coef(lmSpikeIn)["(Intercept)"]
      abline(lm(spikeInDrugControls[,i] ~ spikein$logCellNumber))
   }
   legend("bottom",inset=c(0,-.55),colnames(spikeInDrugControls),cex=0.8)
   dev.off()
}
print("Interpolation Data")
print(head(interpolationInfo))

# Combine DMSO Wells for interpolation
dmsoWells <- gsub("Well\\w+$","",dmsoSamples)
dmsoWells <- unique(dmsoWells)
print("DMSO wells")
print(dmsoWells)

# Now repeat for DMSO samples
for (dmso in dmsoWells){
   print(paste("DMSO:", dmso))
    dmsoFull <-paste(dmso,"Well",sep="")
   spikeInDrug <- spikeInData[,grepl(dmsoFull,colnames(spikeInData))]
   spikeInDrug <- spikeInDrug[order(rownames(spikeInDrug)),]
   spikeInDrugControls <- spikeInDrug[grepl("^sgNT",rownames(spikeInDrug)),]
   spikeInDrugControls <- spikeInDrugControls[6:10,]
   if(identical(dmso,dmsoWells[1])){
      print("Rownames of selected Controls:")
      print(toString(rownames(spikeInDrugControls)))
      print("List of Controls from spike in file:")
      print(toString(controlGuides))
      print("If the two lists above don't match, this will lead to incorrect assumptions.")
   }
   pdf(paste("plots/",controlCell,".sgNT.",dmso,".plot.pdf",sep=""))
   par(mar=c(12.1,4.1,2.1,6))
   matplot(spikein$logCellNumber,spikeInDrugControls,ylab="log2(readCounts)",xlab="log2(cellNumber)",pch=1,main=paste(controlCell,"sgNT",dmso))
   nc<-ncol(spikeInDrugControls)
   print(spikeInDrugControls)
   for (i in 1:nc){
      sample <- colnames(spikeInDrugControls)[i]
      print(sample)
      lmSpikeIn <- lm(spikeInDrugControls[,i] ~ spikein$logCellNumber)
      interpolationInfo[sample,"slope"] <-coef(lmSpikeIn)["spikein$logCellNumber"]
      interpolationInfo[sample,"intercept"] <- coef(lmSpikeIn)["(Intercept)"]
      abline(lm(spikeInDrugControls[,i] ~ spikein$logCellNumber))
   }
   legend("bottom",inset=c(0,-.55),colnames(spikeInDrugControls),cex=0.8)
   dev.off()
}
print("Interpolation Data")
print(head(interpolationInfo))

sgRNAset <- rownames(logNumbers)
sgRNAset <- unique(gsub("^.*\\.","",sgRNAset))
sgRNAset <- sort(sgRNAset)
print("sgRNAs:")
print(sgRNAset)

# Remove the sgNTs from the sgRNAset
sgRNAset <- sgRNAset[grep("^sgNT",sgRNAset,invert=TRUE)]
sgRNAset <- c("sgNT0-4",sgRNAset)
print("Removed controls that are no longer relevant")
print("sgRNAs:")
print(sgRNAset)

numdrugs <- length(drugs)
numsgrnas <- length(sgRNAset)
numcelllines <- length(cellLines)

# Calculate the estimated cell numbers for each barcode pair.
cellNumbers <- data.frame()
for (well in colnames(logNumbers)){
    for (sgRNA in sgRNAset){
    	for (cellLine in cellLines){
	   sgRNAandCell <- paste(cellLine,sgRNA,sep=".")
   	   cellNumbers[sgRNAandCell,well] <- (logNumbers[sgRNAandCell,well]-interpolationInfo[well,"intercept"])/interpolationInfo[well,"slope"]
    	}
    }
}

# Print intermediary numbers, for debugging.
if (intermediary == 1){
   for (cellLine in cellLines){
      for (drug in drugs){
      	 subsetCellNumbers<-cellNumbers[grep(cellLine,rownames(cellNumbers)),grep(drug,colnames(cellNumbers))]
      	 write.table(subsetCellNumbers,paste("intermediary/",cellLine,".",drug,".log2cellNumbers.txt",sep=""))
      }
   }
}

# Remove DMSO from the drugs list
drugs <- drugs[grep("DMSO",drugs,invert=TRUE)]
numdrugs<-length(drugs)
print("nonDMSO drugs:")
print(drugs)

# Calculate median DMSO cell numbers
DMSOmedian<-vector()
for (sgRNA in sgRNAset){
   for (cellLine in cellLines){
       sgRNAandCell <- paste(cellLine,sgRNA,sep=".")
       DMSOsamples <- unlist(cellNumbers[sgRNAandCell,grep("DMSO",colnames(cellNumbers))])
       DMSOsamples <- DMSOsamples[is.finite(DMSOsamples)]
       DMSOmedian[sgRNAandCell]<-median(DMSOsamples)
   }
}
if (intermediary == 1){  
    write.table(DMSOmedian,"intermediary/DMSOmedians.txt")
}
# Subtract median DMSO cell numbers from drug cell numbers for each drug, dose,
# rep, sgRNA combination.
cellSurvival<-data.frame()
for (sgRNA in sgRNAset){
   for (cellLine in cellLines){
       sgRNAandCell <- paste(cellLine,sgRNA,sep=".")
       for (drug in drugs){
           for (dose in colnames(doses)){
              for (rep in replabels){
		          sample <- paste(drug,dose,rep,sep="_")
		          cellSurvival[sgRNAandCell,sample] <- ((2^cellNumbers[sgRNAandCell,sample])/(2^DMSOmedian[sgRNAandCell]))*100
              }
           }    
       }
   }
}

for (cellLine in cellLines){
   for (drug in drugs){
      subsetCellNumbers<-cellSurvival[grep(cellLine,rownames(cellSurvival)),grep(drug,colnames(cellSurvival))]
      write.table(subsetCellNumbers,paste(cellLine,drug,"cellSurvival.txt",sep="."))
   }
}

meanCellSurvival<-data.frame()

# Calculate the mean value for the two replicates for each drug dose for
# each sgRNA.
for (sgRNA in sgRNAset){
  for (cellLine in cellLines){
     sgRNAandCell <- paste(cellLine,sgRNA,sep=".")
     for (drug in drugs){
        for (dose in colnames(doses)){
            rep1<-paste(drug,dose,"Rep1",sep="_")
            rep2<-paste(drug,dose,"Rep2",sep="_")
            newCol<-paste(drug,dose,sep="_")
            meanCellSurvival[sgRNAandCell,newCol]<-mean(c(cellSurvival[sgRNAandCell,rep1],cellSurvival[sgRNAandCell,rep2]))
            #print(paste(sgRNAandCell,rep1,cellSurvival[sgRNAandCell,rep1],rep2,cellSurvival[sgRNAandCell,rep2],meanCellSurvival[sgRNAandCell,newCol],sep="  "))
        }
     }
  }
}

print("Plotting heatmaps.")
my_palette<- colorRampPalette(c("black","green"))(n=199)
col_breaks = c(seq(0,95,length=150),seq(96,150,length=25),seq(156,200,length=25))
for(drug in drugs){
   samples <- as.matrix(data.frame(meanCellSurvival[order(rownames(meanCellSurvival)),grepl(drug,colnames(meanCellSurvival))]))
   drugDose <- doses[drug,] + 6
   drugDose <- 10^(drugDose)
   colnames(samples)<-drugDose
   samples<-samples[,order(colnames(samples))]
   colnames(samples)<-paste(colnames(samples),"uM")
   baseline <- vector(mode="integer", length=numsgrnas)
   baseline <- rep(100,numsgrnas)
   samples <- cbind(baseline,samples)
   write.table(samples,paste(drug,"cellSurvival.heatmap.txt",sep="."))
   pdf(paste("plots/",drug,".cellSurvival.heatmap.pdf",sep=""),width=8.5,height=11)
   heatmap.2(samples,Colv=NA,Rowv=NA,margins=c(6,12),main=paste(drug,"Mean Cell Survival",sep="\n"),col=my_palette,trace="none",dendrogram="none",scale="none",density.info="none",breaks=col_breaks)
   dev.off()
}

# Now for DMSO
samples <- as.matrix(data.frame(justNumbers[order(rownames(justNumbers)),grepl("DMSO",colnames(justNumbers))]))
colnames(samples)<-gsub("\\_\\w+$","",colnames(samples))
print(head(samples))
write.table(samples,paste("DMSO","ReadCounts.heatmap.txt",sep="."))
pdf(paste("plots/DMSO","ReadCounts.heatmap.pdf",sep="."),width=14,height=9)
heatmap.2(samples,Colv=NA,Rowv=NA,margins=c(7,7),main=paste("DMSO","Read Counts",sep="\n"),trace="none",dendrogram="none")
dev.off()

samples <- as.matrix(data.frame(cellNumbers[order(rownames(cellNumbers)),grepl("DMSO",colnames(cellNumbers))]))
colnames(samples)<-gsub("\\_\\w+$","",colnames(samples))
print(head(samples))
write.table(samples,paste("DMSO","CellNumbers.heatmap.txt",sep="."))
pdf(paste("plots/DMSO",".cellNumbers.heatmap.pdf",sep=""),width=14,height=9)
heatmap.2(samples,Colv=NA,Rowv=NA,margins=c(7,7),main=paste("DMSO","Cell Numbers",sep="\n"),trace="none",dendrogram="none")
dev.off()

print("Plotting dose response curve.")
 color_palette<-rainbow(numsgrnas)
# Trying to use the R package drc...
 for(cellLine in cellLines){
    for(drug in drugs){
        print(paste("Plotting data for",cellLine,"and",drug))
        # Pull out the samples for this cell line and drug
        samples <- cellSurvival[grepl(cellLine,rownames(cellSurvival)),grepl(drug,colnames(cellSurvival))]
        # Strip the drug and cell line names out of the row and column labels
        colnames(samples)<-gsub(drug,"",colnames(samples))
        colnames(samples)<-gsub("^_","",colnames(samples))
        rownames(samples)<-gsub(cellLine,"",rownames(samples))
        rownames(samples)<-gsub("^.","",rownames(samples))
        # Look up the doses for this particular drug
        drugDose<-doses[drug,]
        drugDose<-10^drugDose
        # Define a new column for drug concentration (two values per dose for the two reps)
        Concentration<-vector()
        for(i in 1:(numdoses*numreps)){
        	  j <- floor((i+1)/numreps)
 	  Concentration[i] <- drugDose[j]
        }
        # Transpose the sample table
        samples<-t(samples)
        # Add Concentration as the first column in the sample table.
        samples <- data.frame(cbind(Concentration,samples))
        print(samples)
        # Sort from lowest to highest concentration.
#        samples <- samples[order(samples$Concentration),]
        # Write that table to a file.
       # write.table(samples,paste("intermediary/",cellLine,".",drug,".cellSurvival.txt",sep=""))
        # Create a pdf.
 #       pdf(paste("plots/",cellLine,".",drug,".cellSurvival.pdf",sep=""))
 #       colnames(samples)<-gsub("sgNT0.4","sgNT0-4",colnames(samples))
        # Draw a line for each sgRNA, each with a different color in the rainbow.
 #       plotStarted <- 0
#        for (i in 1:numsgrnas){
           # Pick out Concentration column and column for one sgRNA
# 	   subset <- cbind(samples$Concentration,samples[,sgRNAset[i]])
# 	   # Label this data subset and make it a data frame
# 	   rownames(subset)<-rownames(samples)
# 	   colnames(subset)<-c("Concentration","value")
# 	   subset<-data.frame(subset)
 	   # Fit a model. (Use llogistic in drc)
    # Constrain parameter top to 100.
    # model = log(inhibitor) vs response (three parameters)
    # x is the log (dose) and curve will have a standard slope (hill slope = -1.0)
# 	   samples.m1 <- drm(value~Concentration, data = subset, fct = LL.4(),na.action=na.omit,control=drmc(errorm=FALSE))
# 	   print("Fit")
# 	   print(samples.m1$fit)
 	   #write(unlist(samples.m1$fit),paste("intermediary/",cellLine,".",drug,".",sgRNAset[i],"model.txt",sep=""))
 	   # Plot the model.  If the plot hasn't been started, start a new plot. Otherwise add to the plot.
# 	   if (!is.null(samples.m1$fit)){
# 	     if (plotStarted == 0){		
# 	   	plot(samples.m1,type="bars",col=color_palette[i],ylim=c(0,150))
# 		plotStarted <- 1
# 	     } else {   
 	      	#plot(samples.m1,type="bars",col=color_palette[i],ylim=c(0,150),add=TRUE)
# 	     }
# 	    }
# 	}
# 	# This attempt at a legend doesn't work.
# 	#legend("topright",inset=c(-.3,0),sgRNAset,cex=0.8,xpd=TRUE,fill=color_palette)
 	# Close the pdf, end the plot
 #	dev.off()
 #    }
     }
 }
 print("Finished")
