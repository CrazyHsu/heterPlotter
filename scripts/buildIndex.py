#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
File name: buildIndex.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2018-10-15 13:52:33
Last modified: 2018-10-15 13:52:35
'''

import argparse
from commonFunc import *

#================parse arguments============
parser = argparse.ArgumentParser()
parser.add_argument("-bf", "--bigfile", default="example.hapmap", type=str,
                    help="This is the hapmap file you want to build index.")
parser.add_argument("-idx", "--index", default=1000, type=int,
                    help="The number of SNPs you want to specified in each bin. Default: 1000.")
parser.add_argument("-ft", "--fileType", default=1, type=int,
                    help="The number of the file type you want to build index.'1' refer to hapmap file,"
                         " and '2' refer to gwas file.")
parser.add_argument("-o", "--out", default="example.idx", type=str,
                    help="This is the index file you specified.")
#===========================================

args = parser.parse_args()
assert validateFile(args.bigfile)
bigFile = args.bigfile
fileType = args.fileType
out = args.out

buildIndex(bigFile, out, fileType=fileType, lineOffset=args.index)
