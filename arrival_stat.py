#!/usr/bin/env python
# -*- coding: utf-8 -*-



"""
Created on Wed Mar  6 12:41:43 2019

@author: savytska
"""



import sys
import argparse
import os
import pandas as pd
import numpy as np
import math
import statistics as st
import matplotlib.pyplot as plt
from matplotlib import axes
from matplotlib.pyplot import figure
from subprocess import check_output















#############################################################################################################

## Get a value from stadard file, from a specific line

## Obtains the whole line from the file
     
def get_one_line(filepath, line_number):
	return check_output([
		"sed",
		"-n",
		"%sp" % line_number,
		filepath
	])

## Obtains the number from the file

def get_n(filepath,line_number):
	readline=get_one_line(filepath,line_number)
	readline=str(readline)
	readline=readline.replace("b'SN\\treads mapped:\\t","")
	readline=readline.replace("\\n'","")
	return readline

###############################################################################################################




## Obtain the overrepresented in sequencing data regions


## Arrival-like Statistics is calculated for each position in the reference mitogenome (can be closely related species, if absent). For this the per base coverage file is generated with bedtools genecoverage,
## based on the reads mapped to reference (using minimap2, e.g.). 
## To calculate mean coverage for short reads the following formula is used: c=nreads*lreads/lmito, where nreads - n of mapped to ref reads, lreads- length of short reads(hard-coded to 151 in script), lmito-length of mitogenome
## To calculate mean coverage for long reads the following formula is used: c=sum of lengths of all mapped reads/length of mitogenome
## While there are no expected tandem repeats in mammalian mitogenome, that thus could be collapsed by the assembly, there are nevertheless NUMTs(nuclears 

def arrival_st(coverage, cmean):
	coverage.columns=["ref","base","coverage"]
	rep_frag=coverage[coverage["coverage"]>=cmean//math.log(2)]
	if len(rep_frag)==0:
		print("Regions with negative values for Arrival Statistics are not found!\n \
                       Maximal per base coverage: %s \n \
                       Minimal per base coverage: %s \n \
                       Mean coverage: %s"  %(max(coverage["coverage"]),min(coverage["coverage"]),cmean))
		file = open('%s_%s.log' %(output_id,rtype), "a")
		file.write("Regions with negative values for Arrival Statistics are not found!\n \
                       Maximal per base coverage: %s \n \
                       Minimal per base coverage: %s \n \
                       Mean coverage: %s"  %(max(coverage["coverage"]),min(coverage["coverage"]),cmean))
		file.close()
	else:
		print("Regions with negative Arrival Statistics found! Warning! Possible overrepresentation of fragmented or NUMT-originated sequences in the sample!")
		file = open('%s_%s.log' %(output_id,rtype), "a")
		file.write("Regions with negative Arrival Statistics found! Warning! Possible overrepresentation of fragmented or NUMT-originated sequences in the sample!\n \
                       Maximal per base coverage: %s \n \
                       Minimal per base coverage: %s \n \
                       Mean coverage: %s"  %(max(coverage["coverage"]),min(coverage["coverage"]),cmean))
		file.close()
		rep_frag.to_csv("%s.log" %(output_id), sep='\t', mode='a')



                     
## Plotting graphs


def plot_coverage(coverage, cmean,output_id,rtype):
	coverage.columns=["ref","base","coverage"]
	x = coverage["base"]
	y = coverage["coverage"]
	print(x.head())
	print(y.head())

	colors = coverage["coverage"]
	fig, ax = plt.subplots()
	fig.set_figheight(15)
	fig.set_figwidth(25)
	im = ax.scatter(x, y, c=colors, 
            cmap='viridis', label="coverage/base")
	cb=fig.colorbar(im, ax=ax, cmap='viridis')  # show color scale
	plt.xlabel('Reference Mitochondrion, base position', fontsize=34, labelpad=15)
	plt.ylabel('Coverage, read per base', fontsize=34, labelpad=15)
	plt.tick_params(labelsize=25)
	cb.ax.tick_params(labelsize=25)
	ax.axhline(y=cmean, c='red', ls='--', label="mean coverage")
	ax.axhline(y=cmean//math.log(2), c='blue', label="mean coverage/ln2")
	ax.legend()
	print(cmean//math.log(2))
	print(cmean)
	plt.savefig('%s_%s.png' %(output_id,rtype)) 













## Parsing input arguments


#coverage.bed -r ["short","long"] -as align.stat -if index.fai -o output_folder
parser=argparse.ArgumentParser(description="Small program for calculating A stats, based on the newly sequenced short or long reads mapped to the reference sequence. \n \
				For both short and long read sequencing data requires file with per base coverage information, obtained from bedtools genomecov .\n \
				Example of command to obtain per base coverage:        'bedtools genomecov -ibam reads_mapped_to_ref.bam -d > output_cov_pb.bed' \n \
				For short read data the mapping statistics file is also needed, obtained by calling 'samtools stats' on your initial .bam file \n \
				For long read data you need an index file instead, summerizing length of each mapped long read per read id. \n \
				Such index file can be obtained when .bam with mapped reads only is converted back to .fq file and further indexed with samtools fqidx \n \
				The custom name for output files can also be provided with flag -o,--output")
parser.add_argument("coverage", help = "Direct to the file, containing per base coverage information, must be string", type = str)
parser.add_argument("rtype", help = "Read type; 'short', 'long' ", type = str)
parser.add_argument("-as","--alstats", help = "For short read data: alignment statistics file; needed to extract total number of mapped reads. Needed for short reads only", type=str)
parser.add_argument("-if","--indexfile", help = "Fasq index file, if it is long reads that are mapped. Needed for coverage calculation", type=str)
parser.add_argument("-o","--output", help = "Output folder name, must be a string", type=str, default = "NW_alignment.fasta")
args = parser.parse_args()


coverage = args.coverage
rtype = args.rtype
alstats = args.alstats
ifile = args.indexfile
output_id = args.output
######## MAIN #######



def main():

	print("Initializing assessment of the coverage statistics of the mitochondrial assembly, based on the reference")
	
	file = open('%s_%s.log' %(output_id,rtype), "a")
	file.write("Main input file: %s \n \
		   Sequencing read type analyzed: %s \n " %(coverage,rtype))
	file.close()
	if str(rtype)=="short":
		cov_sr=pd.read_csv(str(coverage), sep='\t',header=None)
		nmapped=int(get_n(alstats,14))
		cmean=151*nmapped//len(cov_sr)
		arrival_st(cov_sr, cmean)
		plot_coverage(cov_sr, cmean,output_id,rtype)
	elif str(rtype)=="long":
		faidx_lr=pd.read_csv(str(ifile),sep='\t',header=None)
		cov_lr=pd.read_csv(str(coverage), sep='\t',header=None)
		print(cov_lr.head())
		faidx_lr.columns=["read","len","1","2","3","4"]
		cmean=int(sum(faidx_lr["len"]))//len(cov_lr)
		arrival_st(cov_lr, cmean)
		plot_coverage(cov_lr, cmean,output_id,rtype)
	else:
		sys.exit("Type of the reads was provided incorrectly! Please, provide the suitable arguments in your command! Valid options: short, long")

	
if __name__ == "__main__":
	    main()
