#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
# author Bo Yan

'''
# Logs:

Created on Jan 15, 2020

Modifications on Feg 12, 2020:
    Add DensityPlot function to plot correlation between two tsscount groups
Modifications on April 20, 2020:
    Add absolute path file for output, input
'''

'''
Based on python 2.7

Include multiple handling for comparison between different Cappable-seq tsscount groups, comparison between Enrich vs Control.
Since using temp file based on time, could be called from different threads

Usage:

Compare enrich and control tsscount files, filter TSS positions with EnrichRatio and TPM cutoff
$python /mnt/home/ettwiller/yan/Bo_script/CompareTSS.py filter --input enrich.tsscount.gtf control.tsscount.gtf --output .tsscount.gtf --tpmcutoff --ratiocutoff
    Ordering of input files: enrich.tsscount.gtf control.tsscount.gtf

Compare enrich and control tsscount files, plot x-axis TPM vs y-axis EnrichRatio for TSS positions (above TPM cutoff) in enrich.
$python CompareTSS.py plotEnrich --input enrich.tsscount.gtf control.tsscount.gtf --output .png --tpmcutoff --ratiocutoff
    Ordering of input files: enrich.tsscount.gtf control.tsscount.gtf

Extract the nio/tpm of TSS positions in multiple input files
$python CompareTSS.py extract --input input1.tsscount.gtf input2 input3 --output --verbose optional

Compute and plot the pearson correlation between TSS in input1 and input2, e.g. correlation and scatter plot between replicates
$python CompareTSS.py correlation --input input1.tsscount.gtf input2.tsscount.gtf --output png --tpmcutoff
'''

try:
    import os, sys
    from subprocess import check_call
    import time
    import re
    import argparse
    import math
    from scipy.stats import gaussian_kde
    import numpy as np
    import random
    import matplotlib.pyplot as plt
    import pandas as pd
    plt.switch_backend('agg') # need this for running on server
except:
    print "module error."
    quit()

def CreateTemp():
    '''
    create a temp file based on the time
    '''
    localtime = time.asctime(time.localtime())
    Prefix = ''.join(localtime.split()[-2].split(':'))
    return 'comparison.temp' + Prefix 

def FindTPM_1(input_gtf):
    '''
    Use to find the value of TPM_1 in the input_gtf
    input e.g.:
    NC_000913.3     M9.control.tss  TSS     148     148     9.76438538076   +       1-coordination  nio=2;TPM=9.76438538076;
    
    return 9.76438538076/2, which is a float
    '''
    with open(input_gtf) as f:
        line = f.readline().split('\t')
        nio = int(re.findall('nio=(\d.*?);', line[-1])[0]) 
        TPM = float(re.findall('TPM=(\d.*?);', line[-1])[0]) 
    return TPM/nio


class Pair:
   
    def __init__(self, filelist):
        self.inputlist = filelist 

        self.comparison = None 
        self.TPM_1_control = None
        self.TPM_1_input = None
    

    def compareEnvsCon(self, output_file):
        '''
        Compare Enrich vs Control, enrich should be the first file in self.inputlist.

        Only output the positions with TPM>0 in enrich, and calculate the EnrichRatio,
        for positions having 0 read in control, control_TPM=0

        EnrichRatio=(TPM_Enrich)/(TPM_Control),
        and for the positions having 0 read in control (nio_control=0;TPM_control=0), EnrichRatio=(TPM_Enrich)/(TPM_1_control);
        here I do not round the EnrichRatio, since low fold e.g. 0.003 is rouned to 0.0, which affects the log10(Ratio) transformation.
        In FilterTssGTF.py output, the EnrichRatio is rounded to 2 digits.

        output self.comparison:
        <chr,TSS,strand><enrich_TPM><control_TPM><enrich_nio><control_nio><EnrichRatio>
        '''

        if len(self.inputlist)!=2:
            print "Error: CompareTwo function can take only two input files."
            quit()

        self.TPM_1_control = FindTPM_1(self.inputlist[1])

        with open(self.inputlist[1]) as f1:
            key = [(line.strip().split('\t')[0], line.strip().split('\t')[3], line.strip().split('\t')[6]) for line in f1] # key: (chr, TSS, strand)
            with open(self.inputlist[1]) as f2:
                tpm = (re.findall('TPM=(\d.*?);', line.strip())[0] for line in f2)
                dic_control_tpm = dict(zip(key, tpm))
                f2.seek(0)
                nio = (re.findall('nio=(\d.*?);', line.strip())[0] for line in f2)
                dic_control_nio = dict(zip(key, nio))

        output = open(output_file, 'w')
        count_enrich = 0
        with open(self.inputlist[0]) as f:
            for line in f:
                key = (line.strip().split('\t')[0], line.strip().split('\t')[3], line.strip().split('\t')[6])
                temp = [';'.join(key)]
                TPM_Enrich = re.findall('TPM=(\d.*?);', line.strip())[0]
                nio_Enrich = re.findall('nio=(\d.*?);', line.strip())[0]
                count_enrich += int(nio_Enrich)
                temp.extend([TPM_Enrich, dic_control_tpm.get(key, '0'), nio_Enrich, dic_control_nio.get(key, '0')])
                EnrichRatio = float(TPM_Enrich)/float(dic_control_tpm.get(key, self.TPM_1_control))
                temp.append(str(EnrichRatio))
                print>>output, '\t'.join(temp)
        output.close()
        print "The total number of reads in Enrich is: {}".format(count_enrich)

        self.comparison = output_file

        return


    def calculatePositions(self, TPM_cutoff, EnrichRatio_cutoff):
        '''
        calculate the number of positions and reads with TPM and EnrichRatio cutoff for comparison between enrich and control.
        '''
        if self.comparison:
            data = pd.read_csv(self.comparison, header=None, sep='\t', names = ['Enrich', 'Control', 'EnrichNio', 'ControlNio', 'EnrichRatio'])
            tpm = data.loc[data['Enrich']>=TPM_cutoff, :]
            ratio = data.loc[(data['Enrich']>=TPM_cutoff) & (data['EnrichRatio']>=EnrichRatio_cutoff), :]
            print "Total number of TSS in Enrich: {}.".format(data.shape[0])
            print "Total number of TSS in Enrich with TPM>={}: {}".format(TPM_cutoff, tpm.shape[0])
            print "Total number of TSS in Enrich with TPM>={} and EnrichRatio>={}: {}".format(TPM_cutoff, EnrichRatio_cutoff, ratio.shape[0])
            print "Total TPM of TSS in Enrich: {}.".format(sum(data['Enrich']))
            print "Total TPM of TSS in Enrich with TPM>={}: {}".format(TPM_cutoff, sum(tpm['Enrich']))
            print "Total TPM of TSS in Enrich with TPM>={} and EnrichRatio>={}: {}".format(TPM_cutoff, EnrichRatio_cutoff, sum(ratio['Enrich']))
        else:
            print "Need to run object.compareEnvsCon first."
        return 


    def outputEnvsCon(self, output_file, TPM_cutoff, EnrichRatio_cutoff):
        '''
        Output positions in Enrich with TPM and EnrichRatio >=cutoff.
        The same as FilterTssGTF.py, the only difference is Not round EnrichRatio here.

        Here EnrichRatio: Original division, No round
        self.comparison
        <chr,TSS,strand><enrich_TPM><control_TPM><enrich_nio><control_nio><EnrichRatio>
        '''
        if self.comparison:
            print "Compare enrich and control."
            print "Save positions with TPM>={} and EnrichRatio>={} in output.".format(TPM_cutoff, EnrichRatio_cutoff)

            output = open(output_file, 'w')
            with open(self.comparison) as f:
                for line in f:
                    line = line.strip().split('\t')
                    if float(line[5]) >= EnrichRatio_cutoff and float(line[1]) >= TPM_cutoff:
                        temp = [line[0].split(';')[0], '.', 'tssgtf', line[0].split(';')[1], line[0].split(';')[1], '.', line[0].split(';')[2], '1-coordination']
                        attri = ['nio='+line[3],'TPM='+line[1],'Ratio='+line[5], 'nio_control='+line[4], 'TPM_control='+line[2],''] # I add '' here to have the ';' at the end of '\t'.join(temp)
                        temp.append(';'.join(attri))
                        print>>output, '\t'.join(temp)
            output.close()
        else:
            print "Need to run object.compareEnvsCon first."
        return 


    def jitter(self, list, var):
        '''
        Used to add a value for all the numbers in list for jitter
        value is randomly chosed from normal distribution of (mean=number, variance=var)
        
        var controls the jitter amount
        '''
        ls_jitter = []
        for item in list:
            ls_jitter.append(random.choice(np.random.normal(item, var, size=50)))
        return ls_jitter

    def PlotEnrichRatio(self, png, TPM_cutoff):
        '''
        Plot y-axis EnrichRatio vs x-axis EnrichTPM
        input e.g.
        <chr,TSS,strand><enrich_TPM><control_TPM><EnrichRatio>
        
        TPM_cutoff: only plot the positions with TPM_Enrich>=cutoff
        '''
        if self.comparison:
            print "Compare enrich and control."
            print "Plot y-axis EnrichRatio vs x-axis EnrichTPM."
            data = pd.read_csv(self.comparison, header=None, sep='\t', names = ['Enrich', 'Control', 'EnrichNio', 'ControlNio', 'EnrichRatio'])
            df = data.loc[data['Enrich']>=TPM_cutoff, :]

            x = self.jitter([np.log10(item) for item in df['Enrich']], 0.06) 
            y = self.jitter([np.log10(item) for item in df['EnrichRatio']], 0.06)

            plt.rcParams['font.size'] = 8 

            fig, ax = plt.subplots(dpi=300)
            pcm = ax.scatter(x, y, s=1.5, edgecolor='', alpha=0.6) 
            
            plt.ylabel('log10(EnrichRatio)')
            plt.xlabel('log10(Enrich TPM)')
            ax.set_aspect(1 / ax.get_data_ratio())
            plt.savefig(png, transparent=True)
        else:
            print "Need to run object.compareEnvsCon first."
        return

    def compareGroups(self, output, nio=None):
        '''
        extract TPM of TSS for each input files listed in the input_file (a list).
        Use tpm=0 for positions that do not exist in any input file
        
        output e.g.
        <chr;TSS 1-cooraltion;strand><TPM in input1><TPM in input2><TPM in input3>.
        '''
        
        number_file = len(self.inputlist)
        dic_tpm = {}
        dic_nio = {}
        i = 0
        ls_TPM_1 = [] 

        for input in self.inputlist:
            ls_TPM_1.append(FindTPM_1(input))

            with open(input) as f:
                for line in f:
                    key = (line.strip().split('\t')[0], line.strip().split('\t')[3], line.strip().split('\t')[6])
                    if not dic_tpm.has_key(key): 
                        dic_tpm[key] = ['0']*number_file 
                    dic_tpm[key][i] = re.findall('TPM=(\d.*?);', line.strip().split('\t')[-1])[0]

                    if not dic_nio.has_key(key):
                        dic_nio[key] = ['0']*number_file
                    dic_nio[key][i] = re.findall('nio=(\d.*?);', line.strip().split('\t')[-1])[0]

            i +=1

        with open(output,'w') as f:
            if nio:
                for item in dic_nio:
                    ls = [';'.join([item[0], item[1], item[2]])]
                    ls.extend(dic_nio[item])
                    print>>f, '\t'.join(ls)
            else:
                for item in dic_tpm:
                    ls = [';'.join([item[0], item[1], item[2]])]
                    ls.extend(dic_tpm[item])
                    print>>f, '\t'.join(ls)
        
        self.TPM_1_input = ls_TPM_1[:]
        return ls_TPM_1
    
    def DensityPlot(self, input, cutoff, png):
        '''
        Generate a DensityDot Plot of TPM_input1 vs TPM_input2
        
        input is the output generated by compareGroups, e.g.
        <chr;TSS 1-cooraltion;strand><TPM in input1><TPM in input2>
        
        cutoff: only plot the positions with TPM_input1>=cutoff or TPM_input2>=cutoff
        
        png: name of output png file
        '''
        name = ['rep1', 'rep2']
        data = pd.read_csv(input, header=None, sep='\t', names = [name[0], name[1]])
        
        print "Correlation between two:"
        print np.corrcoef(data[name[0]], data[name[1]])

        df = data.loc[(data[name[0]]>=cutoff) | (data[name[1]]>=cutoff)]

        x = np.asarray([math.log10(self.TPM_1_input[0]) if item ==0 else math.log10(item) for item in df[name[0]]])
        y = np.asarray([math.log10(self.TPM_1_input[1]) if item ==0 else math.log10(item) for item in df[name[1]]])

        xscale = int(math.floor(np.ndarray.max(x)))
        yscale = int(math.floor(np.ndarray.max(y)))

        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

        plt.rcParams['font.size'] = 8 

        fig, ax = plt.subplots(dpi=300)
        pcm = ax.scatter(x, y, c=z, s=2, cmap='coolwarm', edgecolor='', alpha=0.8)

        xaxis, yaxis = [-1], [-1] 
        xaxis.extend(range(0, xscale+1))
        yaxis.extend(range(0, yscale+1))
        xaxis_ticks, yaxis_ticks  = ['0.1'], ['0.1']
        xaxis_ticks.extend([str(int(math.pow(10, item))) for item in range(0, xscale+1)])
        yaxis_ticks.extend([str(int(math.pow(10, item))) for item in range(0, yscale+1)])
        
        plt.xticks(xaxis, xaxis_ticks)
        plt.yticks(yaxis, yaxis_ticks)
        
        plt.xlabel('TPM of {}'.format(name[0]), fontsize=12)
        plt.ylabel('TPM of {}'.format(name[1]), fontsize=12)

        plt.colorbar(pcm, aspect=10, shrink = 0.3)       
        ax.set_aspect(1 / ax.get_data_ratio()) 
        plt.savefig(png, transparent=True) 


##-------Parser
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help', dest = 'mode')
    
    parser_a = subparsers.add_parser('filter', help='compare enrich and control, filter based on TPM and EnrichRatio cutoff')
    parser_a.add_argument('--input', nargs = '+', help = 'tsscount files', dest='input_file')
    parser_a.add_argument('--output', help='output tsscount.gtf containing enrich and control', dest='output_file')
    parser_a.add_argument('--tpmcutoff', help='tpm cutoff for enrich', dest='TPM_cutoff', type=float, default=1)
    parser_a.add_argument('--ratiocutoff', help='EnrichRatio cutoff for enrich', dest='EnrichRatio_cutoff', type=float, default=1)

    parser_b = subparsers.add_parser('plotEnrich', help='compare enrich and control, plot the EnrichRatio and calculate the TSS and Reads')
    parser_b.add_argument('--input', nargs = '+', help = 'tsscount files', dest='input_file')
    parser_b.add_argument('--output', help='png file for EnrichRatio', dest='output_file')
    parser_b.add_argument('--tpmcutoff', help='tpm cutoff for plotting enrich and calculation', dest='TPM_cutoff', type=float, default=1)
    parser_b.add_argument('--ratiocutoff', help='EnrichRatio cutoff for calculation', dest='EnrichRatio_cutoff', type=float, default=1)

    parser_c = subparsers.add_parser('extract', help='extract TSS from different groups')
    parser_c.add_argument('--input', nargs = '+', help = 'tsscount files', dest='input_file')
    parser_c.add_argument('--output', help='files saving TPM or nio for different groups', dest='output_file')
    parser_c.add_argument('--verbose', help='add to output nio instead of tpm', action='store_true')

    parser_d = subparsers.add_parser('correlation', help='extract TSS from repeats and plot correlation')
    parser_d.add_argument('--input', nargs = '+', help = 'tsscount files', dest='input_file')
    parser_d.add_argument('--output', help='png file of correlation', dest='output_file')
    parser_d.add_argument('--tpmcutoff', help='tpm cutoff for plotting two repeats', dest='TPM_cutoff', type=float, default=1)

    args = parser.parse_args()

    Comparelist = Pair([os.path.abspath(item) for item in args.input_file])
    output = os.path.abspath(args.output_file)

    if args.mode =='filter':
        tempfile = CreateTemp()
        Comparelist.compareEnvsCon(tempfile)
        Comparelist.outputEnvsCon(output, args.TPM_cutoff, args.EnrichRatio_cutoff)
        os.remove(tempfile)

    elif args.mode == 'plotEnrich':
        tempfile = CreateTemp()
        Comparelist.compareEnvsCon(tempfile)
        Comparelist.PlotEnrichRatio(output, args.TPM_cutoff)
        Comparelist.calculatePositions(args.TPM_cutoff, args.EnrichRatio_cutoff)
        os.remove(tempfile)
    
    elif args.mode == 'extract':
        Comparelist.compareGroups(output, args.verbose)
    
    elif args.mode == 'correlation':
        temp_output = CreateTemp()
        Comparelist.compareGroups(temp_output)
        Comparelist.DensityPlot(temp_output, args.TPM_cutoff, output)
        os.remove(temp_output)


