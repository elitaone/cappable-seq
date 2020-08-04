#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
# author Bo Yan

'''
# Logs:

Created on Sept 4, 2018

Made modification on Feb 4,2019: 
    change to generate temporary dir based on the current time instead of tempSuffix (therefore remove --sufix), 
    so multiple threads can run this script in parallel without conflicts. 

Made modification on Aug 15,2019:
    remove the mention of bedtools path since does not need bedtools
'''

'''
Based on python 2.7

This is used to count the TPM or nio of each TSS.

Usage:
$python CountTssGTF.py --input .gtf --output tsscount.gtf --cutoff (default 0)

Logic:
For each position in the genome, the algorithm counts the number of Tss tags in each orientation (nio);
output coordination is 1-based.
Since the Tss position is based on the bed file which is converted from the bam file, 
whether using soft clip or not during the mapping is important for the Tss caculation.

--input: a gtf file having the read information, which is generated from bam or bed from Read1 (small RNA library prep) using bam2gtf.py.
         e.g.
         chr1    B_En.primary.bam        gtf     1407172 1407247 .       -       1-coordination  NS500355:NS500355:HLVYYBGXC:2:13105:8846:2837;255

--output: a sorted standard gtf tss file correponding to TSS genomic position. The coordinates are also in 1-based system.
        e.g.
        chr, source, feature, TSS, TSS, TPM, strand, coordination, attribute (nio=number of reads; TPM)
        chr1    B_En.gtf        tssgtf  1407247 1407247 5.18105185715   -       1-coordination  nio=1;TPM=5.18105185715;
        chr1    B_En.gtf        tssgtf  4508552 4508552 5.18105185715   +       1-coordination  nio=1;TPM=5.18105185715;
        Here soure is the input file name.
        For the attribute, I add ';' at the end, which is used for re.findall search in the other functions.

--cutoff: a float used to filter the positions whose TPM are below the defined cutoff (default 0, no filtering) in the output

Note:
(1) This script could be used for counting genome with multiple chrs.
(2) The tmp file has unique name based on computer time, so this script could be called from different threads.
(3) Result is the same as generated by Laurence's bam2firstbasegtf.pl;
    nio: number of reads starting at this TSS positions;
    TPM = nio/(total number of entries in the bed file converted by bamtobed)*1,000,000
(4) output has the precise TSS information that could be used for getfasta directly to have the nucleotide at the TSS position:
    $bedtools getfasta -s -fi /mnt/home/ettwiller/yan/reference/GRC38_hg20/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed .gtf -fo temp.fasta
    >1:15425884-15425885(+) = IGV 15425885
    G
(5) The output gtf is sorted based on the chr and start (sort -f -k1,1 -k4,4n).
'''

try:
    import os
    import re
    import argparse
    from subprocess import check_call, check_output
    import time
except:
    print "Module Error!"
    quit()

##-------functions

def firstbase(input_gtf):
    '''
    Use to create a gtf file containing only the TSS position;
    input split.1, e.g.
    chr1    B_En.primary.bam        gtf     1407172 1407247 .       -       1-coordination  NS500355:NS500355:HLVYYBGXC:2:13105:8846:2837;255
    output split.1.onebase, e.g.
    chr1	B_En.primary.bam	gtf	1407247	1407247	.	-	1-coordination	NS500355:NS500355:HLVYYBGXC:2:13105:8846:2837;255
    '''
    output = open(input_gtf+'.onebase','w')
    with open(input_gtf) as f:
        for line in f:
            line = line.strip().split('\t')
            if line[6] =="+":
                line[4] = line[3]
            else:
                line[3] = line[4]
            print>>output, "\t".join(line)
    output.close()
        
        
def main_countreads(input_gtf, output_gtf, cutoff):
    
    print "============="
    localtime = time.asctime(time.localtime()) 
    
    print "Start at:", localtime
    print "Counting the number of TSS tags for input file {}.".format(os.path.basename(input_gtf))
    
    feature = 'tssgtf'
    source = os.path.basename(input_gtf)
    
    dir = os.getcwd() 
    
    suffix = localtime.split()[-2].split(':')
    suffix.append(localtime.split()[-1]) 
    suffix= ''.join(suffix) 
    print "Creating temporary dir {}".format('temp' + suffix)
    
    try:
        os.mkdir('temp'+suffix)
    except:
        print "Temporary dir exists!"
        quit()

    command = 'cp {} {}'.format(input_gtf, 'temp' + suffix +'/count.temp.gtf')
    check_call(command, shell = True)
    os.chdir('temp'+suffix)
    command = 'wc -l count.temp.gtf'
    totalnumber = int(check_output(command, shell = True).split()[0])
        
    print "There are %d number of reads used for caculating Tss tags." % totalnumber
    command = "awk -F \'\\t' \'{print > $1}\' count.temp.gtf"   
    check_call(command, shell = True)
    
    for item in os.listdir(os.getcwd()):
        if item != 'count.temp.gtf':
            new_name = 'split.' + re.sub('[\|\&]', '_', item) 
            os.rename(item, new_name)
    
    for file in os.listdir(os.getcwd()): 
        if 'split' in os.path.basename(file):
            firstbase(file) 
            command = 'cut -f1-5,7 {}.onebase | sort -k4,4n -k6 | uniq -c > {}.uniqc'.format(file, file)
            check_call(command, shell = True)
    
    output=open(os.path.basename(output_gtf),'w')
    totalTpm = 0
    totalTss = 0
    for file in os.listdir(os.getcwd()):
        if 'uniqc' in os.path.basename(file):
            with open(file) as f:
                for line in f:
                    totalTss +=1
                    list = [line.strip().split()[1], source, feature]
                    list.extend(line.strip().split()[4:6])
                    TPM = int(line.strip().split()[0])*1000000.0/totalnumber
                    list.append(str(TPM))
                    list.append(line.strip().split()[6])
                    list.append('1-coordination')
                    temp = ";".join(['nio=' + line.strip().split()[0], 'TPM='+str(TPM),''])
                    list.append(temp)
                    if TPM>=cutoff:
                        print>>output, '\t'.join(list)
                        totalTpm +=1
                         
    print "There are %d unique TSS tags in the output file." % totalTss
    print "Among these unique TSS tags, there are %d having TPM above %f." % (totalTpm, cutoff)
    output.close()
    
    if os.path.dirname(output_gtf): 
        command = 'sort -k1,1 -k4,4n ' + os.path.basename(output_gtf) + '>' + os.path.abspath(output_gtf)
        dir_output = os.path.dirname(output_gtf)
        if dir_output == '.': 
            dir_output = dir
    else:
        command = 'sort -k1,1 -k4,4n ' + os.path.basename(output_gtf) + '> ../' + os.path.basename(output_gtf)
        dir_output = dir

    check_call(command, shell = True)
    
    print "The output file is {}.".format(os.path.basename(output_gtf))
    print "The output path is {}.".format(dir_output)

    os.chdir(dir)
    command = 'rm -r ' + ' temp' + suffix
    check_call(command, shell = True)

    localtime = time.asctime(time.localtime())
    print "End at :", localtime
    print "============="
    return 1

    
##-------Parser
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input file name', dest='input_file')
    parser.add_argument('-o', '--output', help='output file name', dest='output_file')
    parser.add_argument('-c', '--cutoff', help='cutoff of TPM', type=float, dest='cutoff', default=0)

    args = parser.parse_args()
    
    main_countreads(args.input_file, args.output_file, args.cutoff)
   
    
 
