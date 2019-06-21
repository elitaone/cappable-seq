#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
# author Bo Yan

'''
# Logs:

Created on Sept 3, 2018

Made modifications on Jan 31, 2019: 
    create temporary files based on the current time, so multiple threads can run this script in parallel without conflicts.
'''

'''
Based on python 2.7

This is used to convert the bam/bed file (0-coordination) into a gtf file.

Usage:
$python bam2gtf.py --input bam  --output .gtf --type bam 
OR
$python bam2gtf.py --input bed  --output .gtf --type bed

Logic:
To convert bed (bamtobed) to gtf, I use bed-column2+1 -> gtf-start, bed-column3 -> gtf-end, 
and combine the information in bed-column4(ID)-column5(MAPQ) (if not empty) to the attribute column in gtf (using ';' as seperator).

--input is a bed or sorted and indexed bam file.
bed e.g.
NC_000913.3     0       87      M01193:350:000000000-AAM45:1:1108:18550:9187    44      +

--ouput is a gtf file, containing the read ID and MAPQ
gtf e.g.
NC_000913.3     Replicate1_enriched.bam gtf     1       87      .       +       1-coordination  M01193:350:000000000-AAM45:1:1108:18550:9187;44
Note:
(1) here 44 is from the fifth column of bed file, which is the MAPQ from sam file.

--type: the type of input file, could be bam or bed, default is bam

Note:
(1) Need bedtools in $PATH or change the bedtools path manually in the parser.
(2) bam to bed conversion is performed using bedtools bamtobed, 
    so Flag 4 entry is not in the converted bed file;
    do not do any Flag filter of reads in bam file,
    therefore need to perform Flag selection (e.g. -F 256) before running this conversion step.
'''

try:
    import argparse
    import sys
    import os
    from subprocess import check_call
    import time
except:
    print "Module Error!"
    quit()


#//----------
def bamtobed(input_bam, path):
    '''
    convert bam to bed using bedtools;
    the output bed file is suffix.count.temp.bed;
    note I do not do any filter of the reads in bam file, 
    therefore perform filter (e.g. Flag) before this step if necessary.
    '''
    localtime = time.asctime(time.localtime())
    
    # create temporary file name, which based on local time, so it will not be the same for different threads, e.g. temp1129362019
    suffix = localtime.split()[-2].split(':')
    suffix.append(localtime.split()[-1]) 
    suffix= ''.join(suffix) # 1533232019
    
    if os.path.isfile('bamtobed.{}.sh'.format(suffix)) or os.path.isfile('{}.count.temp.bed'.format(suffix)):
        print "Tempory file exists! Quit!"
        quit()
    
    with open('bamtobed.{}.sh'.format(suffix), 'w') as script:
        print>>script, '#!/bin/sh'
        print>>script, '{}bedtools bamtobed -i {} > {}'.format(path, input_bam, '{}.count.temp.bed'.format(suffix))
    
    try:
        check_call('sh bamtobed.{}.sh'.format(suffix), shell=True)
    except:
        print "bamtobed Error!"
        os.remove('bamtobed.{}.sh'.format(suffix))
        quit()
        
    os.remove('bamtobed.{}.sh'.format(suffix))
    return '{}.count.temp.bed'.format(suffix) # the name of bed file


def bedtogtf(input_file, output_file, source):

    if os.path.dirname(output_file): # e.g. --output /mnt/home/ettwiller/yan/Ira_CappableSeq/analysisfile/Replicate1_enrich.tss
        dir_output = os.path.dirname(output_file)
        if dir_output == '.':
            dir_output = os.getcwd()
    else: # e.g. --output Replicate1_enrich.tss
        dir_output = os.getcwd()
    
    print "=================="
    print "Start at:", time.asctime(time.localtime())
    print "Convert input file into gtf file."
    print "The output is saved at {}.".format(dir_output)
    
    feature = 'gtf'

    output=open(output_file,'w')
    with open(input_file) as f:
        for line in f:
            line = line.strip().split('\t')
            # bed to gtf
            start = str(int(line[1])+1)
            end = line[2]
            strand = line[5]
            if line[4] != '.':
                attribute = line[3]+';'+line[4]
            else:
                attribute =line[3]

            temp = [line[0], source, feature, start, end, '.', strand, '1-coordination', attribute]
            print>>output, '\t'.join(temp)
    output.close()
    print "Done."
    print "=================="
    return 1
        

##-------Parser
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input file name', dest='input_file')
    parser.add_argument('-o', '--output', help='output file name', dest='output_file')
    parser.add_argument('-t', '--type', help='input file type', type=str, dest='type', default='bam')

    args = parser.parse_args()

    source=os.path.basename(args.input_file)
    
    if args.type == 'bam':
        # bedtools path
        path = '' # bedtools calling from $PATH; define bedtools path if it is not in the $PATH
        
        bedfilename = bamtobed(args.input_file, path)
        
        # source is the input, feature is gtf 
        bedtogtf(bedfilename, args.output_file, source)
        os.remove(bedfilename)
        
    elif args.type == 'bed':
        bedtogtf(args.input_file, args.output_file, source)

