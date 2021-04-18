# QMAP-Seq
QMAP-Seq code for converting fastq files into relative cell numbers, plotting dose-response curves, calculating AUCs, and generating chemical-genetic interaction networks.  <a href="https://codeocean.com/capsule/3022355/tree/v1"/>
* Code Ocean capsule</a>
* <a href="https://www.nature.com/articles/s41467-020-19553-8First">Original publication</a> 

First, fastq files are converted into relative cell numbers. 

<img src="https://github.com/mendillolab/QMaPP-Seq/blob/master/qmap_flowchart_1.png?raw=true"/>

Once relative cell numbers are calculated for each combination, the data can be aggregated into one concatenated table (concatenate.R). The concatenated dataframe can then be transformed (tidy.R) and used to identify hits (test.R). 

<img src="https://github.com/mendillolab/QMaPP-Seq/blob/master/qmap_flowchart_2.jpg?raw=true"/>
