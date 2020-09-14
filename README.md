# QMAP-Seq
This repository contains QMAP-Seq code for converting fastq files into relative cell numbers, plotting dose-response curves, calculating AUCs, and generating chemical-genetic interaction networks. 
<center>
<img src="https://github.com/mendillolab/QMaPP-Seq/blob/master/qmap_flowchart_1.png?raw=true"/>
</center>
Once relative cell numbers are calculated for each combination, the data can be aggregated by concatenating each table. Once the data has been concatenated, it can be transformed (tidy.R) and used to identify hits (test.R). 

<img src="https://github.com/mendillolab/QMaPP-Seq/blob/master/qmap_flowchart_2.jpg?raw=true"/>
