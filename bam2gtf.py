#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
# author Bo Yan

'''
# Logs:

Created on Sept 3, 2018

Made modifications on Jan 31, 2019: 
    create temporary files based on the current time, so multiple threads can run this script in parallel without conflicts.
    
Made modifications on June 24, 2019:
    add important information about bamtobed conversion in the Logic part, also add the printing in bamtobed function.
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
    '''
    localtime = time.asctime(time.localtime())
    
    suffix = localtime.split()[-2].split(':')
    suffix.append(localtime.split()[-1]) 
    suffix= ''.join(suffix)
    
    if os.path.isfile('bamtobed.{}.sh'.format(suffix)) or os.path.isfile('{}.count.temp.bed'.format(suffix)):
        print "Tempory file exists! Quit!"
        quit()
    
    with open('bamtobed.{}.sh'.format(suffix), 'w') as script:
        print>>script, '#!/bin/sh'
        print>>script, '{}bedtools bamtobed -i {} > {}'.format(path, input_bam, '{}.count.temp.bed'.format(suffix))
    
    try:
        print "bamtobed! Make sure the bam file contains only the primary reads, or one entry for each read."
        check_call('sh bamtobed.{}.sh'.format(suffix), shell=True)
        os.remove('bamtobed.{}.sh'.format(suffix))
    except:
        print "bamtobed Error!"
        os.remove('bamtobed.{}.sh'.format(suffix))
        quit()
        
    return '{}.count.temp.bed'.format(suffix)


def bedtogtf(input_file, output_file, source):

    if os.path.dirname(output_file): # e.g. --output /mnt/home/ettwiller/yan/Ira_CappableSeq/TSS/Replicate1_enrich.tsscount.gtf
        dir_output = os.path.dirname(output_file)
        if dir_output == '.':
            dir_output = os.getcwd()
    else: # e.g. --output Replicate1_enrich.tsscount.gtf
        dir_output = os.getcwd()
    
    
    print "Convert bed file into gtf file."
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
    return 1
        

##-------Parser
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='input file name', dest='input_file')
    parser.add_argument('-o', '--output', help='output file name', dest='output_file')
    parser.add_argument('-t', '--type', help='input file type', type=str, dest='type', default='bam')

    args = parser.parse_args()

    source=os.path.basename(args.input_file)
    
    print "=================="
    print "Start at:", time.asctime(time.localtime())
    
    if args.type == 'bam':

        path = '' # bedtools PATH /home/yan/bin/
        
        bedfilename = bamtobed(args.input_file, path)
        
        bedtogtf(bedfilename, args.output_file, source)
        os.remove(bedfilename)
        
    elif args.type == 'bed':
        bedtogtf(args.input_file, args.output_file, source)
    
    print "Done."
    print "=================="
