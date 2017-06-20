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


def caculate_dict(simu_dict, vcf_dict):
    Simu_number = len(simu_dict)
    vcf_number = len(vcf_dict)
    com_pos = set(vcf_dict.keys()) & set(simu_dict.keys())
    PP = len(set(vcf_dict.keys()) & set(simu_dict.keys()))
    PPR = float(PP) / float(vcf_number)
    N_position = 0
    for pos in com_pos:
        # print vcf_dict[pos], simu_dict[pos]
        # if len(simu_dict[pos][1]) == 2:
        simurecord = (simu_dict[pos][0],
                      [x for x in simu_dict[pos][1] if x != simu_dict[pos][0]][0],
                      simu_dict[pos][2],
                      simu_dict[pos][3])
        # print vcf_dict[pos], simurecord, simu_dict[pos]
        # sys.exit(0)
        # print vcf_dict[pos] , simurecord
        if vcf_dict[pos] == simurecord:
            N_position += 1
    N_Accuracy = N_position
    Accuracy_Rate = float(N_Accuracy) / float(Simu_number)
    Positive_rate = float(N_position) / float(vcf_number)
    N_FP = PP - N_position + len(set(vcf_dict.keys())) - PP
    FPR = float(N_FP) / float(vcf_number)
    N_FN = Simu_number - PP
    FNR = float(Simu_number - PP) / float(Simu_number)
    return [
        Simu_number,
        vcf_number,
        N_Accuracy,
        Accuracy_Rate,
        PP,
        PPR,
        N_position,
        Positive_rate,
        N_FP,
        FPR,
        N_FN,
        FNR]


def caculate_info(index, length, depth):
    vcf_file = './simu_%s_%sX_record_vcf_picklefile.txt' % (length, depth)
    simu_file = './simu_%s_%sX_record_simu_picklefile.txt' % (length, depth)
    output_file = './simu_%s_%sX_record_statics.txt' % (length, depth)
    pickle_tuple = pickle.load(open(vcf_file, 'rb'))
    vcf_dict, repeatregion_vcf_dict, non_repeatregion_vcf_dict = pickle_tuple[
        0], pickle_tuple[1], pickle_tuple[2]
    pickle_tuple = pickle.load(open(simu_file, 'rb'))
    simu_dict, repeatregion_simu_dict, non_repeatregion_simu_dict = pickle_tuple[
        0], pickle_tuple[1], pickle_tuple[2]
    output_handle = csv.writer(open(output_file, "w"), delimiter="\t")
    output_handle.writerow(["Sample",
                            "Length",
                            "Depth",
                            "#Simu",
                            "#Find",
                            "Accuracy",
                            "Accuracy_Rate",
                            "#PP",
                            "PPR",
                            "#Positive",
                            "Positive_Rate",
                            "#FP",
                            "FPR",
                            "#FN",
                            "FNR"])
    output_handle.writerow([index, length, depth] +
                           caculate_dict(simu_dict, vcf_dict))
    output_handle.writerow([index, length, depth] +
                           caculate_dict(repeatregion_simu_dict, repeatregion_vcf_dict))
    output_handle.writerow([index,
                            length,
                            depth] + caculate_dict(non_repeatregion_simu_dict,
                                                   non_repeatregion_vcf_dict))


def main():
    length, depth, index = prj_name.split("_")[1], prj_name.split("_")[
        2][0:-1], prj_name.split("_")[-1]

    caculate_info(index, length, depth)


if __name__ == '__main__':
    prj_folder = os.getcwd()
    prj_tree = ProjectFolders(prj_folder)
    prj_name = fullpath2last_folder(prj_tree.home)
    main()
