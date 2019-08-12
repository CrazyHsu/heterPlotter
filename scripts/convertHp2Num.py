#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
File name: convertHp2Num.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2018-10-18 21:08:23
Last modified: 2018-10-18 21:08:24
'''

import argparse
from commonFunc import *

#================parse arguments============
parser = argparse.ArgumentParser()
parser.add_argument("-hp", "--hapmap", default="example.hapmap", type=str,
                    help="This is the config file which contains the absolute path of hapmap files or"
                         "a part section of hapmap file.")
#===========================================

args = parser.parse_args()
assert validateFile(args.hapmap)
hapmapList = getFileList(args.hapmap)
curDir = os.getcwd()
genoTypeDict = {}
for i in range(len(hapmapList)):
    hapmap = hapmapList[i]
    assert hapmap
    hapmapBase = os.path.basename(hapmap)
    hapmapNum = os.path.join(curDir, hapmapBase + ".num")
    hpFileTmpPect = os.path.join(curDir, hapmapBase + ".pect")
    if i == 0:
        genoTypeDict = convertHapmap2Num_new(hapmap, hapmapNum, hpFileTmpPect, genoTypeDict)
    else:
        convertHapmap2Num_new(hapmap, hapmapNum, hpFileTmpPect, genoTypeDict)

