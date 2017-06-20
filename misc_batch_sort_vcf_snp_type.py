#!/usr/bin/env python
# encoding: utf-8
"""
miss_trans_nucltide_on_hg19.py

Created by Mingchen on 2014-12-11.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile:  original fasta file
orgnism:  human,mouse or rabbit
referecce_file:  the file contain variable region start and end info eg: human_get_vdj.txt

note: already trans the reversed reads
"""
import sys
import os
import csv
import re
import glob
import copy
import subprocess
import time
import multiprocessing
import traceback
import tempfile
import numpy as np
import matplotlib.pyplot as plt
from stackedBarGraph import StackedBarGrapher
import Bio.Alphabet
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, Process, Manager
from bsub import bsub
from mytools import *
from misc_prepare_pbs import *
from misc_get_trimmed_region import *
from itertools import groupby
try:
    import cPickle as pickle
except ImportError:
    import pickle
from collections import Counter


def get_simu_dict(length, depth):
    repeat_region = csv.reader(
        open(
            "/zzh_gpfs02/zhangyanfang/repeat_region_analysis/all_repeat_region_interval.txt",
            "rU"),
        delimiter="\t")
    repeat_region_dict = {}
    for index, line in enumerate(repeat_region):
        repeat_region_dict.setdefault(line[0], []).append((line[1], line[2]))

    simu_record = csv.reader(
        open(
            "./simu_%s_%sX_record.txt" %
            (length, depth), "rU"), delimiter="\t")

    simu_dict, repeatregion_simu_dict, non_repeatregion_simu_dict = {}, {}, {}
    print IUB_CODE
    for index, line in enumerate(simu_record):
        if index % 10000 == 0:
            print "Processed 10000 lines in simu %s %s." % (length, depth)
        # print line
        chromsome = line[0].split("_")[0]
        position = int(line[1]) + int(line[0].split("_")[1])
        ref_nucle = line[2]
        snv_nucle = IUB_CODE[line[3]]
        haplo_type = line[4]
        in_repeat_region = 0
        # print chromsome, position, ref_nucle, snv_nucle, haplo_type,
        # in_repeat_region
        try:
            for region in repeat_region_dict[chromsome]:
                if int(region[0]) <= position - 1 <= int(region[1]):
                    in_repeat_region = 1
                    # print chromsome, position, ref_nucle, snv_nucle, haplo_type, in_repeat_region
                    # sys.exit(0)
                    break
        except KeyError:
            pass

        # print ref_nucle, snv_nucle, haplo_type
        simu_dict[(chromsome, position)] = (
            ref_nucle, snv_nucle, haplo_type, in_repeat_region)
        if in_repeat_region == 1:
            repeatregion_simu_dict[(chromsome, position)] = (
                ref_nucle, snv_nucle, haplo_type, in_repeat_region)
        else:
            non_repeatregion_simu_dict[(chromsome, position)] = (
                ref_nucle, snv_nucle, haplo_type, in_repeat_region)
        # print simu_dict
        # if index == 10:
        #	sys.exit(0)

    # How many position be detected
    # print len(set(vcf_dict.keys())), len(set(simu_dict.keys())), len(set(vcf_dict.keys()) & set(simu_dict.keys()))
    # print len(set(repeatregion_vcf_dict.keys())), len(set(repeatregion_simu_dict.keys())), len(set(repeatregion_vcf_dict.keys()) & set(repeatregion_simu_dict.keys()))
    # print len(set(non_repeatregion_vcf_dict.keys())),
    # len(set(non_repeatregion_simu_dict.keys())),
    # len(set(non_repeatregion_vcf_dict.keys()) &
    # set(non_repeatregion_simu_dict.keys()))
    pickle_file = './simu_%s_%sX_record_simu_picklefile.txt' % (length, depth)
    pickle_file_handle = open(pickle_file, 'wb')
    dump_tuple = (
        simu_dict,
        repeatregion_simu_dict,
        non_repeatregion_simu_dict)
    pickle.dump(dump_tuple, pickle_file_handle)
    pickle_file_handle.close()


def get_vcf_dict(length, depth):
    print "processing length:%s, depth:%s" % (length, depth)
    repeat_region = csv.reader(
        open(
            "/zzh_gpfs02/zhangyanfang/repeat_region_analysis/all_repeat_region_interval.txt",
            "rU"),
        delimiter="\t")
    repeat_region_dict = {}
    for index, line in enumerate(repeat_region):
        repeat_region_dict.setdefault(line[0], []).append((line[1], line[2]))

    vcf_file = "./6-haplorcaller-wholegenome-member/simu_%s_%sX.filter.uniq.concordant.merge.sort.dedup.AddGroup.snp.haplo.g.vcf.haplercaller.vcf" % (
        length, depth)
    handle = csv.reader(open(vcf_file, "rU"), delimiter="\t")
    vcf_dict, repeatregion_vcf_dict, non_repeatregion_vcf_dict, haplo_type_set = {}, {}, {}, []
    haplo_type_r = []
    haplodict = {"1/1": "-", "0/1": "+", "1/0": "+", "0/0": "-"}
    for index, line in enumerate(handle):
        if index % 10000 == 0:
            print "Processed 10000 lines in simu %s %s." % (length, depth)
        if "chr" in line[0] and len(line) == 10:
            chromsome = line[0]
            position = int(line[1])
            ref_nucle = line[3]
            snv_nucle = line[4]
            haplo_type_r.append(line[9].split(":")[0])
            in_repeat_region = 0
            try:
                for region in repeat_region_dict[chromsome]:
                    if int(region[0]) <= position - 1 <= int(region[1]):
                        in_repeat_region = 1
            except KeyError:
                pass
            try:
                haplo_type = haplodict[line[9].split(":")[0]]
            except KeyError:
                haplo_type = "-"
            # print chromsome, ref_nucle, snv_nucle, haplo_type
            vcf_dict[(chromsome, position)] = (
                ref_nucle, snv_nucle, haplo_type, in_repeat_region)
            if in_repeat_region == 1:
                repeatregion_vcf_dict[(chromsome, position)] = (
                    ref_nucle, snv_nucle, haplo_type, in_repeat_region)
            else:
                non_repeatregion_vcf_dict[(chromsome, position)] = (
                    ref_nucle, snv_nucle, haplo_type, in_repeat_region)
            haplo_type_set.append(haplo_type)
            # sys.exit(0)
    #'''
    # print Counter(haplo_type_set), Counter(haplo_type_r)
    # How many position be detected
    # print len(set(vcf_dict.keys())), len(set(simu_dict.keys())), len(set(vcf_dict.keys()) & set(simu_dict.keys()))
    # print len(set(repeatregion_vcf_dict.keys())), len(set(repeatregion_simu_dict.keys())), len(set(repeatregion_vcf_dict.keys()) & set(repeatregion_simu_dict.keys()))
    # print len(set(non_repeatregion_vcf_dict.keys())),
    # len(set(non_repeatregion_simu_dict.keys())),
    # len(set(non_repeatregion_vcf_dict.keys()) &
    # set(non_repeatregion_simu_dict.keys()))
    pickle_file = './simu_%s_%sX_record_vcf_picklefile.txt' % (length, depth)
    pickle_file_handle = open(pickle_file, 'wb')
    dump_tuple = (vcf_dict, repeatregion_vcf_dict, non_repeatregion_vcf_dict)
    pickle.dump(dump_tuple, pickle_file_handle)
    pickle_file_handle.close()


def main():
    length, depth, index = prj_name.split("_")[1], prj_name.split("_")[
        2][0:-1], prj_name.split("_")[-1]

    get_vcf_dict(length, depth)
    get_simu_dict(length, depth)


    #caculate_info(index, length, depth)
if __name__ == '__main__':
    prj_folder = os.getcwd()
    prj_tree = ProjectFolders(prj_folder)
    prj_name = fullpath2last_folder(prj_tree.home)
    main()
