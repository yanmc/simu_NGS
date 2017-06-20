#!/usr/bin/env python
# encoding: utf-8
"""
3.0.py

Created by Mingchen on 2015-05-15.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.
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
    for file_name in os.listdir(prj_folder):
        print file_name
        # and( "K" in file_name or "L" in file_name):
        if os.path.isdir(file_name):
            os.chdir(file_name)
            print os.getcwd()  # , os.system("which 2.0.py")
            if len(file_name.split("_")) == 4 and int(
                    file_name.split("_")[-1]) <= 20:
                print "Yes!"
                length, depth = file_name.split(
                    "_")[1], file_name.split("_")[2][0:-1]
                print length, depth
                reads_num = int(depth) * 3.27 * 100000
                #os.system("bsub -n 1 -q zzh -e err.0 -o out.0 \"~/simulator/wgsim/wgsim -e 0 -1 %s -2 %s -r 0.001 -R 0 -X 0 -S -1 -N %s /zzh_gpfs02/zhangyanfang/repeat_region_analysis/Long_identical_region_up_down_same_length.fa simu_%s_%sX.read1.fq simu_%s_%sX.read2.fq > simu_%s_%sX_record.txt &\""%(int(length), int(length), reads_num, length, depth, length, depth, length, depth))
                # print "~/simulator/wgsim/wgsim -e 0 -1 %s -2 %s -r 0.001 -R 0 -X 0 -S -1 -N %s /zzh_gpfs02/zhangyanfang/repeat_region_analysis/Long_identical_region_up_down_1000bp.fa simu_%s_%sX.read1.fq simu_%s_%sX.read2.fq > simu_%s_%sX_record.txt & "%(int(length), int(length), reads_num, length, depth, length, depth, length, depth)
                #os.system("cp /zzh_gpfs02/yanmingchen/WGS_simulation/repeat_region_simu/simu_150_100X/*.py ./")
                #os.system("python whole_script.py -i %s/\*.fq -seq wholeseq -qc 20 &"%(os.getcwd()))
                #os.system("rm -rf ./1-trim-file")
                #os.system("rm -rf ./2-bwa-alignment-file")
                #os.system("rm -rf ./3-filter-file")
                #os.system("bsub -n 1 -q zzh -e err.0 -o out.0 \"misc_batch_sort_vcf_snp_type.py\"")
                os.system(
                    "bsub -n 1 -q zzh -e err.0 -o out.0 \"misc_batch_caculate_vcf_accuracy.py\"")
            os.chdir(prj_folder)


if __name__ == '__main__':

    pool_size = multiprocessing.cpu_count()
    # create 1st and 2nd subfolders
    prj_folder = os.getcwd()
    prj_tree = ProjectFolders(prj_folder)
    prj_name = fullpath2last_folder(prj_tree.home)
    start = time.time()

    main()
    end = time.time()
    print prj_name, end - start
    print "Finished"
