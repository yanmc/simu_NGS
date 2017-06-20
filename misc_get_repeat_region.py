#!/usr/bin/env python
# encoding: utf-8
"""
3.0.py 

Created by Mingchen on 2015-05-15.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.
"""

import sys, os, csv, re, glob, copy, subprocess, time, multiprocessing
import traceback, tempfile
import pandas as pd
import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
#from stackedBarGraph import StackedBarGrapher
import Bio.Alphabet 
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, Process, Manager
from bsub import bsub
from mytools import *
#from misc_prepare_pbs import *
#from misc_get_trimmed_region import *
from collections import Counter
try:
    import cPickle as pickle
except ImportError:
    import pickle
import statsmodels.api as sm

def main():
	repeat_region_file = "/zzh_gpfs02/zhangyanfang/repeat_region_analysis/all_repeat_region_interval.txt"
	genome_file = "/zzh_gpfs/data/ref_genome/hg19_genome/hg19.fa"
	gennome_record_handle = SeqIO.parse(open(genome_file, "rU"), "fasta")
	gennome_record = {}
	for record in gennome_record_handle:
		gennome_record[record.id] = record
	
	handle = csv.reader(open(snv_records_file, "rU"), delimiter = "\t")
	handle.next()
	handle.next()
	info, result = [], {}
	for index, line in enumerate(handle):
		if index %1000 == 0:
			print "%s has Done!"%index
		query_chr = line[0]
		query_position = int(line[4]) -1
		subject_chr = line[6]
		subject_position =  int(line[9]) -1
		ref_nucle = line[10]
		alt_nucle = line[11]
		
		

if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	start = time.time()

	main()
	end = time.time()
	print prj_name, end-start
	print "Finished"
