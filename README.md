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
> Step 1: Pre-processing of the raw fastq file to select the reads arising from 3' end of the transcripts. Following are the processing steps:<br/>
  * Getting the reverse complement of reads
  * Trimming the terminal polyA to retain only two A's in the end (or begining)
  * Excluding those reads that contains N (undefined nucleotide) within them
  * Filtering out the sequences that are less than 20 nucleotides in length <br/>
<code> ./3Pseq_iniprocess.py -h </code><br/>
<code>usage: 3Pseq_iniprocess.py [-h] -q FASTQ [-g GENOME] [-c CONFIG] [-o OUTPUT] </code><br/>

Script performs the initial processing <br/>
optional arguments: <br/>
-h, --help -> show this help message and exit <br/>
-q FASTQ, --fastq FASTQ ->  RAW FASTQ file of the reads <br/>
-g GENOME, --genome GENOME -> Path to the complete genome file <br/>
-c CONFIG, --config CONFIG -> Configuration file <br/>
-o OUTPUT, --output OUTPUT -> Output directory to store the results <br/>

- - - -

> Step 2: Aligning the processed raw FASTQ file to the genome. The script performs the following actions:
  * Aligning the reads to the genome using <a href="http://bowtie-bio.sourceforge.net/index.shtml" target="_blank">bowtie</a>.</br>
  * Filtering out the alignment to retain those reads that contain the 'A' mismatch either at the end or the beginning (This would ensure it is the untemplated 'A' that was added during the process of cleavage/polyadenylation.)<br/>
  * Give the bed output of the alignment which will later be used for deriving the 3P-Peaks <br/>
<code> ./alignment_trigger.py -h</code><br/>
<code>usage: alignment_trigger.py [-h] [-q FASTQ] -g GENOME -c CONFIG [-o OUTPUT]</code><br/>

Script perfoms the alignment and filtering <br/>
optional arguments: <br/>
-h, --help -> show this help message and exit <br/>
-q FASTQ, --fastq  -> FASTQ Processed FASTQ file of the reads [output from 3Pseq_iniprocess.py] </br>
-g GENOME, --genome GENOME -> Path to the bowtie genome index file </br>
-c CONFIG, --config CONFIG -> Configuration file </br>
-o OUTPUT, --output OUTPUT -> Output directory to store the results </br>

- - - -

> Step 3: Deriving the 3P-Peaks.<br/>
<code>
Usage: ./detect_peaks \<bedcount_file\> \<readcut-off\>
</code><br/>

Please ensure that input bedcount file for the script is in the following format:

First column - Contig/Chromosome name <br/>
Second column - Strand <br/>
Third column - Alignment begin <br/>
Fourth column - Alignment end <br/>
Fifth column - Read count <br/>
Sixt column - The lenght of the range #Not required <br/>





