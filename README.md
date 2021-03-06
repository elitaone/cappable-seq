## Cappable-seq Analysis
This document explains all the scripts used for analyzing Cappable-seq data. <br>
Scripts are developed and maintained by Bo Yan (New England Biolabs, yan@neb.com). <br>

******************
### **Pipeline for Cappable-seq or ReCappable-seq TSS analysis** <br>
- TSS library is prepared using small RNA library (ligation based) <br>
- trim the adapter using cutadapt for illumina Read1 <br>
- map Read1 using bowtie2 or STAR with soft clipping at the 5'end <br>
- convert bam file to gtf using **bam2gtf.py** <br>
- count the number of reads starting at each TSS (5'end tag) using **CountTssGTP.py** <br>
- compare with the duplicates for correlation using **CompareTSS.py** <br>
- compare enrich to the control library to remove the false positive TSS using **CompareTSS.py** <br>
******************

### **bam2gtf.py** <br>
Convert bam or bed file (0-coordination) into gtf file (1-coordination). <br>
<pre>
python bam2gtf.py --input input.bam --output output.gtf --type bam
python bam2gtf.py --input input.bed --output output.gtf --type bed
</pre>
input: a bed or sorted and indexed bam file. <br>
output: a gtf file. <br>
type: bam or bed, default bam

_Note_: <br>
Run by python2; <br>
Require bedtools pre-installed in PATH, otherwise need to change the path of bedtools; <br>
Remove reads with Flag 4 (umapped reads) in bam file, therefore need to perform other Flag selection before running this conversion step <br>
The attribution (column9, using ';' as delimiter) in gtf output contains both column 4 and column 5 information from bed file, or ReadID and MAPQ score from bam file. e.g. M01193:350:000000000-AAM45:1:1108:18550:9187;44 <br>


### **CountTssGTF.py** <br>
Count the TPM (number of 5'end tag per million of reads) or nio for each TSS. <br>
<br>
nio: number of reads starting at this TSS positions (5'end tag); <br>
TPM = nio/(total number of entries in the bed file converted by bamtobed)*1,000,000
<pre>
<b>python CountTssGTF.py --input input.gtf --output output.tsscount.gtf --cutoff float</b>
</pre>
input: a gtf file having the read information, which is generated from bam or bed from Read1 using bam2gtf.py. <br>
output: a gtf file containing the TPM and nio (tsscount.gtf file) corresponding to TSS genomic position, sorted by chromosome and start. Column5: TPM, Attribution: nio;TPM;. e.g.<br>
<pre>NC_000913.3     Replicate1_enriched.gtf tssgtf  12      12      0.15128805897   +       1-coordination  nio=2;TPM=0.15128805897;</pre> 
cutoff: a float, only the TSS positions with TPM>=cutoff are reported in the output. Default 0 means no cutoff applied. <br>

_Note_: <br>
Run by python2; <br>
Output has the precise TSS position that could be used for getfasta directly to have the nucleotide at the TSS position. <br>


### **CompareTSS.py** <br>
Contain multiple options to compare the TSS between different groups. <br>
<br>
tsscount.gtf: gtf file (tsscount.gtf) containing the TPM and nio that is generated by CountTssGTF.py. <br>
EnrichRatio=(TPM_Enrich)/(TPM_Control) <br>
TPM_1_control: TPM of position having 1 read (nio=1) in control <br>
For the positions having 0 read in control (nio_control=0;TPM_control=0), EnrichRatio=(TPM_Enrich)/(TPM_1_control) <br>

_Note_: <br>
Run by python2; <br> 
Need numpy, scipy.stats, pandas and matplotlib pre-installed. <br>

<div><b>* Option: plotEnrich</b></div>
Plot x-axis Enrich TPM vs y-axis EnrichRatio. Each dot represents a TSS in enrich. A successful enrichment will show two seperate populations of TSS.
<pre><b>
python CompareTSS.py plotEnrich --input enrich.tsscount.gtf control.tsscount.gtf --output output.png --tpmcutoff float --ratiocutoff float
</b></pre>
input: tsscount.gtf file for enrich (first file) and control (second file). <br>
output: png image <br>
tpmcutoff: a float. Only positions in enrich with TPM>=cutoff will be plot. Default 1. <br>
ratiocutoff: a float. Calculate the number of reads and number of TSS positions in enrich with EnrichRatio>=ratiocutoff and TPM>=tpmcutoff. Default 1. <br>
<br>

<div><b>* Option: filter</b></div>
Extract the TSS positions in enrich with EnrichRatio>=ratiocutoff and TPM>=tpmcutoff and save in output. <br>
<pre><b>
python CompareTSS.py filter --input enrich.tsscount.gtf control.tsscount.gtf --output output.tsscount.gtf --tpmcutoff float --ratiocutoff float
</b></pre>
input: tsscount.gtf file for enrich and control. <br>
output: tsscount.gtf file containing TSS positions in enrich above both cutoff. <br>
tpmcutoff: a float. Only positions in enrich with TPM>=tpmcutoff are saved in output. Default 1. <br>
ratiocutoff: a float. Only positions in enrich with EnrichRatio>=ratiocutoff are saved in output. Default 1. <br>

<br>
<div><b>* Option: extract</b></div>
Extract the nio/tpm of all TSS positions from multiple input files. <br>
<pre><b>
python CompareTSS.py extract --input input1.tsscount.gtf input2.tsscount.gtf --output output --verbose optional
</b></pre>
input: tsscount.gtf file for enrich and control. <br>
output: <br>
a gtf file containing the TSS positions for all the input files. Use tpm=0 for positions that do not exist in any input file. <br>
Column1: chr;TSS 1-cooraltion;strand, Column2: TPM in input1, Column3: TPM in input2 ... <br>
verbose: Add this option to report nio instead of TPM in output. <br>

<br>
<div><b>* Option: correlation</b></div>
Compute the pearson correlation between TPM of all the TSS positions in input1 and input2. <br>
Plot (using color to repsent the dot density) the TPM of input1, input2 on x-axis, y-axis. <br>
Only TSS positions with TPM>=tpmcutoff in either input1 or input2 are plot. <br>
<pre><b>
python CompareTSS.py correlation --input input1.tsscount.gtf input2.tsscount.gtf --output output.png --tpmcutoff float
</b></pre>
input1 and input2: tsscount.gtf file for input1 and input2. <br>
output: png image <br>
tpmcutoff: a float. Only TSS positions with TPM>=tpmcutoff in either input1 or input2 are plot. Default 1. <br>

