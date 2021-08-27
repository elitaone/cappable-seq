#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
# author Bo Yan

'''
# Logs:

Created on July 31, 2021
    Combine bam2gtf.py and CountTssGTF.py into one script; add testing bedtools
'''

'''
Based on python2.7
Require bedtools preinstalled

Usage:
$python CountTssFromBam.py --input .bam --output tsscount.gtf --cutoff (default 0) --path (default None) --type (default bam)

--input:
    sorted and indexed bam file contains only primary mapping (samtools -F 256) for Read1.
--output: a sorted standard gtf tss file correponding to TSS genomic position. The coordinates are also in 1-based system.
        e.g.
        chr, source, feature, TSS, TSS, TPM, strand, coordination, attribute (nio=number of reads; TPM)
        chr1    B_En.gtf        tssgtf  1407247 1407247 5.18105185715   -       1-coordination  nio=1;TPM=5.18105185715;
        chr1    B_En.gtf        tssgtf  4508552 4508552 5.18105185715   +       1-coordination  nio=1;TPM=5.18105185715;
        Here source is the input file name.
        nio: number of TSS tags
        TPM: TSS tags per million of mappable reads
--cutoff: a float used to filter the positions whose TPM are below the defined cutoff (default 0, no filtering) in the output
--type: type of input, bam or bed, Default bam
--path: bedtools path, default None means that bedtools is executable in PATH
        e.g. /mnt/home/yan/bin

'''

try:
    import argparse
    import sys
    import os
    from subprocess import check_call, check_output
    import time
    import re
except:
    print "Module Error!"
    quit()

#//----------
def TestCommand(PATH):
    '''
    Test whether bedtools is installed
    quit if any of them can not be found
    '''
    
    with open('testcommand.sh', 'w') as f:
        print>>f, '{}bedtools --version'.format(PATH)
        print>>f, 'status=$?'
        print>>f, 'echo \"returning code: $status\"'
        
    display = check_output('sh testcommand.sh' , shell=True)
    
    code = [item.split(':')[-1].strip() for item in display.split('\n') if item.startswith('returning')] # e.g. [0, 0]
    os.remove('testcommand.sh')
    if code[0] == '0':
        print "Find bedtools."
    else:
        print "Error: Check whether bedtools is installed or PATH is correct."
        quit()

def CheckOutputPath(file):
    if os.path.dirname(file): # e.g. --output /mnt/home/ettwiller/yan/Ira_CappableSeq/TSS/Replicate1_enrich.tsscount.gtf
        dir_output = os.path.abspath(file)
    else: # e.g. --output Replicate1_enrich.tsscount.gtf
        dir_output = os.path.join(os.getcwd(), file)
    return dir_output # the PATH of output_file

def CreatSuffix():
    localtime = time.asctime(time.localtime())
    suffix = localtime.split()[-2].split(':')
    suffix.append(localtime.split()[-1]) 
    suffix= ''.join(suffix)
    return suffix 

#//----------
class Bed:
    def __init__(self, input_bam, output_file, type, PATH):
        self.input_bam = input_bam
        self.output_file = output_file
        self.PATH = PATH # PATH of bedtools
        self.suffix= CreatSuffix()
        
        if type == 'bam':
            tempfile = self.input_file = self.bamtobed()
            self.bedtogtf()
            os.remove(tempfile)
        else:
            self.input_file = self.input_bam
            self.bedtogtf()

    def bamtobed(self):
        '''
        convert bam to bed using bedtools;
        '''
        
        if os.path.isfile('bamtobed.{}.sh'.format(self.suffix)) or os.path.isfile('{}.count.temp.bed'.format(self.suffix)):
            print "Tempory file exists! Quit!"
            quit()
        
        with open('bamtobed.{}.sh'.format(self.suffix), 'w') as script:
            print>>script, '#!/bin/sh'
            print>>script, '{}bedtools bamtobed -i {} > {}'.format(self.PATH, self.input_bam, '{}.count.temp.bed'.format(self.suffix))
        
        try:
            print "bamtobed. Make sure the bam file contains only the primary reads (samtools -F 256)."
            check_call('sh bamtobed.{}.sh'.format(self.suffix), shell=True)
            os.remove('bamtobed.{}.sh'.format(self.suffix))
        except:
            print "bamtobed Error!"
            os.remove('bamtobed.{}.sh'.format(self.suffix))
            quit()
            
        return '{}.count.temp.bed'.format(self.suffix)

    def bedtogtf(self):
        '''Convert bed to gtf'''        
        print "Convert bed file into gtf file."
       
        feature = 'gtf'

        output=open(self.output_file,'w')
        with open(self.input_file) as f:
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
                temp = [line[0], '.', feature, start, end, '.', strand, '1-coordination', attribute]
                print>>output, '\t'.join(temp)
        output.close()
        return self.output_file

#//----------   
class Count:
    def __init__(self, input_gtf, output_file, cutoff, source):
        self.input_gtf = input_gtf
        self.output_file = output_file
        self.source = source # source in the output gtf file
        self.suffix= CreatSuffix()

        self.main_countreads(cutoff)

    def firstbase(self, input_file):
        '''
        Use to create a gtf file containing only the TSS position;
        input split.1, e.g.
        chr1    B_En.primary.bam        gtf     1407172 1407247 .       -       1-coordination  NS500355:NS500355:HLVYYBGXC:2:13105:8846:2837;255
        output split.1.onebase, e.g.
        chr1	B_En.primary.bam	gtf	1407247	1407247	.	-	1-coordination	NS500355:NS500355:HLVYYBGXC:2:13105:8846:2837;255
        '''
        output = open('{}.onebase'.format(input_file),'w')
        with open(input_file) as f:
            for line in f:
                line = line.strip().split('\t')
                if line[6] =="+":
                    line[4] = line[3]
                else:
                    line[3] = line[4]
                print>>output, "\t".join(line)
        output.close()
        return '{}.onebase'.format(input_file)
                
    def main_countreads(self, cutoff):
        
        feature = 'tssgtf'
        dir = os.getcwd() 
      
        try:
            os.mkdir('temp'+suffix)
        except:
            print "Temporary dir exists!"
            quit()

        command = 'cp {} {}'.format(self.input_gtf, 'temp' + suffix +'/count.temp.gtf')
        check_call(command, shell = True)
        os.chdir('temp'+suffix)
        command = 'wc -l count.temp.gtf'
        totalnumber = int(check_output(command, shell = True).split()[0])
            
        print "There are {} number of reads used for caculating TSS tags.".format(totalnumber)
        # Split the file based on the col1 chr: chr1, chr2...
        command = "awk -F \'\\t' \'{print > $1}\' count.temp.gtf"   
        check_call(command, shell = True)
        
        for item in os.listdir(os.getcwd()):
            if item != 'count.temp.gtf':
                new_name = 'split.' + re.sub('[\|\&]', '_', item) # replace '|' in chr e.g. ecoli reference to avoid syntax error
                os.rename(item, new_name)
        
        for file in os.listdir(os.getcwd()): 
            if 'split' in os.path.basename(file):
                self.firstbase(file) 
                command = 'cut -f1-5,7 {}.onebase | sort -k4,4n -k6 | uniq -c > {}.uniqc'.format(file, file)
                check_call(command, shell = True)
        
        output=open('{}.temp'.format(os.path.basename(self.output_file)),'w')
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
        print "Among these unique TSS, there are %d having TPM above %f." % (totalTpm, cutoff)
        output.close()

        # Sort based on chr and TSS
        command = 'sort -f -k1,1 -V -k4,4n {}.temp > {}'.format(os.path.basename(self.output_file), self.output_file)
        check_call(command, shell = True)
        
        print "The output file is {}.".format(self.output_file)

        os.chdir(dir)
        command = 'rm -r ' + ' temp' + suffix
        check_call(command, shell = True)
        return 1

##-------Parser
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input file name', dest='input_bam')
    parser.add_argument('-o', '--output', help='output tsscount file name', dest='output_file')
    parser.add_argument('-t', '--type', help='input file type, bam or bed, default bam', type=str, dest='type', default='bam')
    parser.add_argument('-c', '--cutoff', help='cutoff of TPM', type=float, dest='cutoff', default=0)
    parser.add_argument('-p', '--path', help='Bedtools path, default in PATH', type=str, dest='PATH', default='')

    args = parser.parse_args()

    
    print "=================="
    print "Start at:", time.asctime(time.localtime())
    print "Count the number TSS tags."

    TestCommand(args.PATH)
    # Find abspath for output_file
    output_file = CheckOutputPath(args.output_file)
    
    suffix = CreatSuffix()
    source = os.path.basename(args.input_bam) # source in the output gtf file
    
    Bed(args.input_bam, '{}.{}tempgtf'.format(args.output_file, suffix), args.type, args.PATH)
    # generate gtf file '{}.suffixtemp.gtf'.format(args.output_file) for next step counting

    Count('{}.{}tempgtf'.format(args.output_file, suffix), output_file, args.cutoff, source)
    
    os.remove('{}.{}tempgtf'.format(args.output_file, suffix))

    print "Done."
    print "=================="
