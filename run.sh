#!/usr/bin/env bash
set -ex

pipelineDir=../code/

# Specify location of raw fastq files
rawFastqDir=../data/fastq/

# The rawFastqDir is expected to contain also the Barcodes, DoseList, SpikeInNumbers, etc.  See below.
# maxVar sets the maximum number of mismatches between expected and observed barcodes
maxVar=1
outputDir=../results/

# Set up directories
mkdir $outputDir/concatFastq
mkdir $outputDir/finalCounts

# Iterate through each fastq file in the fastq directory, combining lanes.
for file1 in $rawFastqDir/*L001_R1_001.fastq.gz
do
      echo $file1
      sample=${file1##*/}
      sample=${sample%"_L001_R1_001.fastq.gz"}
      echo $sample
      zcat $rawFastqDir/$sample\_L00*_R1_001.fastq.gz > $outputDir/concatFastq/$sample.CONCAT.fastq
done

# For each concatenated fastq file:
for fastq in $outputDir/concatFastq/*.fastq
do
    # Figure out the sample name from the fastq file name.
    sample=${fastq##*/}
    sample=${sample%".CONCAT.fastq"}
    sample=${sample%"_L001_R1_001.fastq"}
    echo $fastq
    echo $sample

    # Pull out the Barcodes from each read and count the number of times each pair is seen.
    perl $pipelineDir/parseReads2counts.pl $fastq > $outputDir/finalCounts/$sample.rawCounts.txt

    # Read in the raw counts along with the list of real barcodes and label the counts with sgRNA, cell
    (perl $pipelineDir/filterTable.pl $outputDir/finalCounts/$sample.rawCounts.txt $rawFastqDir/Barcodes.txt $maxVar > $outputDir/finalCounts/$sample.labeledCounts.txt) >& $outputDir/finalCounts/$sample.filter.log 
done

# Combine the counts/*labeledCounts.txt files with one column per sample and output one file per cellLine
# Setting intermediary to 1 prints out intermediary data tables.
Rscript $pipelineDir/analyzeData.R --pattern=*labeledCounts.txt --directory=$outputDir/finalCounts --doselist=$rawFastqDir/DoseList.txt --spikein=$rawFastqDir/SpikeInNumbers.txt --intermediary=1

mkdir $outputDir/heatmaps
mkdir $outputDir/cellSurvival
mkdir $outputDir/doseResponseCurves
mv $outputDir/finalCounts/plots $outputDir/plots
mv $outputDir/finalCounts/intermediary $outputDir/intermediary
mv $outputDir/*/*.heatmap.*  $outputDir/heatmaps/
mv $outputDir/*/*.cellSurvival.*  $outputDir/cellSurvival/

