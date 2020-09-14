# QMAP-Seq
This repository contains code for processing data from the QMAP-Seq protocol. First, fastq files are converted into relative cell numbers. 

<img src="https://github.com/mendillolab/QMaPP-Seq/blob/master/qmap_flowchart_1.png?raw=true"/>

Once relative cell numbers are calculated for each combination, the data can be aggregated into one concatenated table. The concatenated dataframe can then be transformed (tidy.R) and used to identify hits (test.R). 

<img src="https://github.com/mendillolab/QMaPP-Seq/blob/master/qmap_flowchart_2.jpg?raw=true"/>
