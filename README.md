## Cappable-seq Analysis
This document explains all the scripts used for analyzing Cappable-seq data. <br>
Scripts are developed and maintained by Bo Yan (New England Biolabs, yan@neb.com). <br>

### **Pipeline for Cappable-seq or ReCappable-seq TSS analysis** <br>
- TSS library is prepared using small RNA library (ligation based) preparation kit <br>
- trim the adapter using cutadapt for illumina Read1 <br>
- map Read1 using bowtie2 or STAR with soft clipping at the 5'end <br>
- count the number of reads starting at each TSS (5'end tag) using **CountTssFromBam.py** <br>
- compare the technical duplicates for correlation using **CompareTSS.py** <br>
- compare enrich to the control library to remove the false positive TSS using **CompareTSS.py** <br>

### Getting Started
---------------

### **CountTssFromBam.py** <br>
Count the TPM (number of 5'end tag per million of reads) or nio for each TSS starting with bam file. <br>
<br>
Note: <br>
Run by python2; <br>
Count all the mappable reads presenting in bam file, therefore select the primary mapping (-F 256, -F 2028) on the STAR bam file if necessary before this step <br>
nio: number of reads starting at this TSS positions (5'end tag); <br>
TPM = nio/(total number of entries in the bed file converted by bamtobed) * 1,000,000 <br>
<pre><b>python CountTssFromBam.py --input input.bam --output output.tsscount.gtf --cutoff float --type bam --path bedtoolsPATH</b> </pre>
input: a sorted and indexed bam file containing the primary reads for TSS counting. or a bed file converteb from the bam file. <br>
output: a gtf file containing the TPM and nio (tsscount.gtf file), sorted by chromosome and start. <br>
<ol>    
Column1, 2, 3, 4, 5, 6, 7, 8: chr, source, feature, TSS, TSS, TPM, strand, coordination <br>
Column9: TPM, Attribution: nio;TPM;. e.g. nio=2;TPM=0.15128805897
</ol>
cutoff: a float, Default 0. Only the TSS positions with TPM>=cutoff are reported in the output. Default 0 means no cutoff applied. <br>
type: type of input, bam or bed, Default bam <br>
path: Default None. bedtools PATH, e.g. ~/bedtools. If bedtools in PATH, leave it as Default. <br>

### **CompareTSS.py** <br>
Contain multiple options to compare the TSS between different groups. <br>
<br>
Note: <br>
Run by python2; <br> 
Need numpy, scipy.stats, pandas and matplotlib module installed. <br>
tsscount.gtf file: file generated by CountTssFromBam.py.
EnrichRatio=(TPM_Enrich)/(TPM_Control) <br>
TPM_1_control: TPM of position having 1 read (nio=1) in control <br>
For the positions having 0 read in control (nio_control=0;TPM_control=0), EnrichRatio=(TPM_Enrich)/(TPM_1_control) <br>

<div><b>* Option: plotEnrich</b></div>
Plot x-axis Enrich TPM vs y-axis EnrichRatio. Each dot represents a TSS in enrich. A successful enrichment will show two seperate populations of TSS.
<pre>
<b>python CompareTSS.py plotEnrich --input enrich.tsscount.gtf control.tsscount.gtf --output output.png --tpmcutoff float --ratiocutoff float</b>
</pre>
input: a tsscount.gtf file for enrich (first file) and control (second file). <br>
output: a png image <br>
tpmcutoff: a float. Default 1. Only positions in enrich with TPM>=cutoff will be plot in png. <br>
ratiocutoff: a float. Default 1. Calculate the number of reads and number of TSS positions in enrich with EnrichRatio>=ratiocutoff and TPM>=tpmcutoff. <br>
<br>

<div><b>* Option: filter</b></div>
Extract the TSS positions in enrich with EnrichRatio>=ratiocutoff and TPM>=tpmcutoff and save in output. <br>
<pre>
<b>python CompareTSS.py filter --input enrich.tsscount.gtf control.tsscount.gtf --output output.tsscount.gtf --tpmcutoff float --ratiocutoff float</b>
</pre>
input: tsscount.gtf file for enrich and control. <br>
output: tsscount.gtf file containing TSS positions in enrich above both cutoff. <br>
tpmcutoff: a float. Default 1. Only positions in enrich with TPM>=tpmcutoff are saved in output. <br>
ratiocutoff: a float. Default 1. Only positions in enrich with EnrichRatio>=ratiocutoff are saved in output. <br>

<br>
<div><b>* Option: extract</b></div>
Extract the nio/tpm of all TSS positions from multiple input files. <br>
<pre>
<b>python CompareTSS.py extract --input input1.tsscount.gtf input2.tsscount.gtf --output output
python CompareTSS.py extract --input input1.tsscount.gtf input2.tsscount.gtf --output output --verbose</b>
</pre>
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
<pre>
<b>python CompareTSS.py correlation --input input1.tsscount.gtf input2.tsscount.gtf --output output.png --tpmcutoff float --Flag AND</b>
</pre>
input: tsscount.gtf file for input1 (x-axis) and input2 (y-axis). <br>
output: a png image <br>
tpmcutoff: a float. Default 1 <br>
Flag: Default AND. Control whether to show the positions only existing in one input in the png image.
<ol>    
Plot TSS positions with TPM>=tpmcutoff in either input1 or input2 with --Flag OR <br>
Plot TSS positions with TPM>=tpmcutoff in both input1 and input2 with --Flag AND <br>
</ol>
