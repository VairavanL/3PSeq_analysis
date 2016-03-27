# Analysis of 3P-Seq data

<p align="center">
<img src="https://github.com/VairavanL/3PSeq_analysis/blob/master/3P_Header.gif"/>
</p>
These are the set of scripts that were used to analyze the 3P-Seq data reported in '***Genome-wide analysis of polyadenylation events in <i>Schmidtea mediterranea</i>***'. These have been tested on sytem with Ubuntu 14.04 LTS Server edition OS. Before you run these scripts, please ensure following dependencies are met:

* Dependencies
  * <a href="http://bowtie-bio.sourceforge.net/index.shtml" target="_blank">Bowtie</a>
  * <a href="http://bedtools.readthedocs.org/en/latest/" target="_blank">Bedtools</a>
  * <a href="http://samtools.sourceforge.net/">Samtools</a>
  * Following python packages:
    * <a href="https://pypi.python.org/pypi/joblib">joblib</a>
    * <a href="http://www.scipy.org/">scipy</a>
    * <a href="http://www.numpy.org/">numpy</a>

###Following is the step-by-step guide to follow for analysis of 3P-Seq data using these scripts :
> Step 1: Pre-processing of the raw fastq file to select the reads arising from 3' end of the transcripts.
 ./3Pseq_iniprocess.py -h <br/>
usage: 3Pseq_iniprocess.py [-h] -q FASTQ [-g GENOME] [-c CONFIG] [-o OUTPUT] <br/>
Script is a wrapper for 3P-Seq pipeline <br/>
optional arguments: <br/>
-h, --help -> show this help message and exit <br/>
-q FASTQ, --fastq FASTQ ->  RAW FASTQ file of the reads <br/>
-g GENOME, --genome GENOME -> Path to the complete genome file <br/>
-c CONFIG, --config CONFIG -> Configuration file <br/>
-o OUTPUT, --output OUTPUT -> Output directory to store the results <br/>
