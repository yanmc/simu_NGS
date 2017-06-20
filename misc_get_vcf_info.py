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
from scipy.stats.stats import pearsonr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import NullFormatter
import math
def get_vcf_dict(length, depth):
	print "processing length:%s, depth:%s"%(length, depth)
	repeat_region = csv.reader(open("/zzh_gpfs02/zhangyanfang/repeat_region_analysis/integrate_by_script/chr_link/new_v3/final_link_new_v6.txt", "rU"), delimiter = "\t")
	repeat_region_dict = {}
	for index, line in enumerate(repeat_region):
		for record in line:
			re_chromesome = "_".join(record.split("_")[:-2])
			re_start = record.split("_")[-2]
			re_end = record.split("_")[-1]
			
			repeat_region_dict.setdefault(re_chromesome, []).append((re_start, re_end))


	vcf_file = "/zzh_gpfs02/yanmingchen/wholegenome/hapler_combine_order/simu_%s_%sX.filter.uniq.concordant.merge.sort.dedup.AddGroup.snp.haplo.g.vcf.haplercaller.vcf"%(length, depth)
	handle = csv.reader(open(vcf_file, "rU"), delimiter = "\t")
	vcf_dict, repeatregion_vcf_dict, non_repeatregion_vcf_dict, haplo_type_set = {},{},{}, []
	haplo_type_r = []
	haplodict = {"1/1": "-",  "0/1": "+", "1/0": "+", "0/0": "-"}
	for index, line in enumerate(handle):
		if index % 100000 == 0:
			print "Processed 10000 lines in simu %s %s."%(length, depth)
		if "chr" in line[0] and len(line) == 10:
			chromsome = line[0]
			position  = int(line[1])
			ref_nucle = line[3]
			snv_nucle = line[4]
			haplo_type_r.append(line[9].split(":")[0])
			in_repeat_region = []
			try:
				for region in repeat_region_dict[chromsome]:
					if int(region[0]) <= position-1 <= int(region[1]):
						#in_repeat_region = "Yes"
						in_repeat_region.append(chromsome + ":" + region[0] + ":" + region[1] + ":" + str(position+1 - int(region[0])) + ":" + str(int(region[1]) - position+1))
			except KeyError:
				pass
			try:
				haplo_type = line[9]
			except KeyError:
				haplo_type = "NA"
			#print chromsome, ref_nucle, snv_nucle, haplo_type
			if in_repeat_region == []:
				vcf_dict[(chromsome, position)] = [chromsome, position, ref_nucle, snv_nucle, haplo_type, "NO"]
			else:
				vcf_dict[(chromsome, position)] = [chromsome, position, ref_nucle, snv_nucle, haplo_type, "|".join(in_repeat_region)]
			haplo_type_set.append(haplo_type)
			#sys.exit(0)
	#'''
	#print Counter(haplo_type_set), Counter(haplo_type_r)
	simu_record_file = "/zzh_gpfs02/yanmingchen/WGS_simulation/simu_%s_%sX_record.txt"%(length, depth)
	simu_handle = csv.reader(open(simu_record_file, "rU"), delimiter = "\t")
	for index, line in enumerate(simu_handle):
		if index % 100000 == 0:
			print "Processed 10000 lines in simu %s %s."%(length, depth)
		#print line
		chromsome = line[0]
		position  = int(line[1])
		ref_nucle = line[2]
		snv_nucle_copy = copy.deepcopy(IUB_CODE[line[3]])
		if len(snv_nucle_copy) == 2:
			snv_nucle_copy.remove(line[2])
			#print snv_nucle_copy, snv_nucle_copy, line[2], type(line[2]), type(snv_nucle_copy)
			snv_nucle = snv_nucle_copy[0]
			
		else:
			snv_nucle = snv_nucle_copy[0]
		haplo_type = line[4]
		try:
			vcf_dict[(chromsome, position)] = vcf_dict[(chromsome, position)] + [ref_nucle, snv_nucle, haplo_type, "Yes", "Yes"]
		except KeyError:
			vcf_dict[(chromsome, position)] = [chromsome, position, "NA", "NA", "NA", "NA"] + [ref_nucle, snv_nucle, haplo_type, "Yes", "No"]
	for (key, value) in vcf_dict.items():
		if len(value) < 9:
			vcf_dict[key] =  vcf_dict[(chromsome, position)] + ["NA", "NA", "NA", "No", "Yes"]
	
	all_one_snv_file = "/zzh_gpfs02/zhangyanfang/repeat_region_analysis/one_snv_analysis/new/all_snv_from_extend_unique.vcf"
	all_one_snv_region_file = "/zzh_gpfs02/zhangyanfang/repeat_region_analysis/one_snv_analysis/new/all_snv_from_extend.vcf"
	all_one_snv_region_handle = csv.reader(open(all_one_snv_region_file, "rU"), delimiter = "\t")
	all_one_snv_region_handle.next()
	all_one_snv_region_handle.next()
	for line in all_one_snv_region_handle:
		chromsome = line[0]
		position  = int(line[4])
		ref_nucle = line[-2]
		snv_nucle = line[-1]
		right_distance = int(line[3]) - int(line[5])
		left_distance = int(line[5])
		distance_record = chromsome + ":" + line[1] + ":" + line[2] + ":" + str(right_distance) + ":" + str(left_distance)
		try:
			if len(vcf_dict[(chromsome, position)]) == 11:
				vcf_dict[(chromsome, position)] = vcf_dict[(chromsome, position)] + [ref_nucle, snv_nucle]
				vcf_dict_dis_record = vcf_dict[(chromsome, position)][5].split("|")
				vcf_dict_dis_record.append(distance_record)
				vcf_dict[(chromsome, position)][5] = "|".join(list(set(vcf_dict_dis_record)))
				
			else:
				vcf_dict[(chromsome, position)][-1] = snv_nucle
				#print vcf_dict[(chromsome, position)][5], type(vcf_dict[(chromsome, position)][5])
				#print vcf_dict[(chromsome, position)][5].split("|") , type(vcf_dict[(chromsome, position)][5].split("|")),distance_record
				distance_record_list = vcf_dict[(chromsome, position)][5].split("|")
				#print distance_record_list, type(distance_record)
				distance_record_list.append(distance_record)
				#print distance_record_list
				#print set(distance_record_list)
				#print list(set(distance_record_list))
				#print "|".join(list(set(distance_record_list)))
				vcf_dict[(chromsome, position)][5] = "|".join(list(set(distance_record_list)))
		except KeyError:
			
			vcf_dict[(chromsome, position)] = [chromsome, position, "NA", "NA", "NA", distance_record,  "NA", "NA", "NA",  "NA", "NA", ref_nucle, snv_nucle]
		
		chromsome = line[6]
		position  = int(line[9])
		ref_nucle = line[-2]
		snv_nucle = line[-1]
		right_distance = int(line[9]) - int(line[7])
		left_distance = int(line[8]) - int(line[9])
		distance_record = chromsome + ":" + line[7] + ":" + line[8] + ":" + str(right_distance) + ":" + str(left_distance)
		try:
			if len(vcf_dict[(chromsome, position)]) == 11:
				vcf_dict[(chromsome, position)] = vcf_dict[(chromsome, position)] + [ref_nucle, snv_nucle]
				vcf_dict_dis_record = vcf_dict[(chromsome, position)][5].split("|")
				vcf_dict_dis_record.append(distance_record)
				vcf_dict[(chromsome, position)][5] = "|".join(list(set(vcf_dict_dis_record)))
			else:
				vcf_dict[(chromsome, position)][-1] = snv_nucle
				vcf_dict_dis_record = vcf_dict[(chromsome, position)][5].split("|")
				vcf_dict_dis_record.append(distance_record)
				vcf_dict[(chromsome, position)][5] = "|".join(list(set(vcf_dict_dis_record)))
		except KeyError:
			
			vcf_dict[(chromsome, position)] = [chromsome, position, "NA", "NA", "NA", distance_record,  "NA", "NA", "NA",  "NA", "NA", ref_nucle, snv_nucle]
		
	#all_one_snv_handle = csv.reader(open(all_one_snv_file, "rU"), delimiter = "\t")
	#for index, line in enumerate(all_one_snv_handle):
	writer = csv.writer(open("%s/%s_%sX_simu_vcf_info.txt"%(prj_folder, length, depth), "w"), delimiter="\t")
	writer.writerow(["chromsome", "position", "ref_nucle", "alt_nucle", "haplo_type", "in_repeat_region", "simu_ref_nucle", "simu_alt_nucle", "simu_haplo_type", "Simu", "VCF_got", "one_snv_ref", "one_snv_simu"])
	writer.writerows(vcf_dict.values())
	
def main():
	#task_pool = Pool(processes = 8)
	index = 0
	for length in [75, 150]:
		for depth in [10, 30, 50]:
			#'''
			index += 1
			#pjobs_ids = task_pool.apply_async(get_vcf_dict, args=( length, depth,))
			get_vcf_dict(length, depth)
	print "Waiting for all subprocesses done..."
	#task_pool.close()
	#task_pool.join()
	
	print 'All subprocesses done.'
	
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
