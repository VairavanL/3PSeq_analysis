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

###Following is the step-by-step guide for analysis of 3P-Seq data using these scripts :
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
Read cut-off -> Minimum number of reads that each nucleotide base should be spanned by.
Please ensure that input bedcount file for the script is in the following format:

First column - Contig/Chromosome name <br/>
Second column - Strand <br/>
Third column - Alignment begin <br/>
Fourth column - Alignment end <br/>
Fifth column - Read count <br/>
Sixt column - The lenght of the range #Not required <br/>

```javascript
Contig22218     -       6342    6390    1       48
Contig22218     -       9626    9675    5451    49
Contig22218     -       9627    9676    55      49
Contig22218     -       9629    9677    42      48
Contig22218     -       9630    9676    13      46
Contig22218     -       9649    9698    1       49
Contig4549      -       41328   41373   20      45
Contig4549      -       41489   41536   71      47
Contig4549      -       41492   41541   205     49
Contig4549      -       41493   41542   140     49
```
- - - -
> Step 4: Detecting the enriched hexameric nucleotide sequence in the 3P-Peak.<br/>
  * The peak sequences can be derived from the output of the previous script. <br/>
  * This scripts will scan through all the peak sequences to search for the conserved hexameric signal using a window-frame approach.<br/>
  * The statistal significance of each hexameric sequence along with its position on the 3P-Peak will be calculated.
<code>
Usage:./hexamer_logic.py \<fastafile\> \<seed\> <br/>
fastafile - The fasta file containing all the 3P-Peak sequences.<br/>
seed - The seed number that will be used for randomization.<br/>
</code>

- - - -

> Step 5: Determining the secondary hexameric PAS signals from the 3P-Peaks <br/>
  * The secondary hexameric PAS signal will be determined from the 3P-Peaks sequences. Fischers test of independence will be used for detecting these signals.<br>
<code>
Usage: ./fischers_test.py \<peak_fasta_file\>  
</code>


