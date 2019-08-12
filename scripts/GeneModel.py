#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
File name: GeneModel.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2018-09-03 22:12:12
Last modified: 2018-09-03 22:12:13
'''

from Features import *

AGPV3_GTF, AGPV3_GFF3, AGPV4_GTF, AGPV4_GFF3 = "", "", "", ""

def getRefGeneModelFile(refGeneModelNo, modelAnnoType, **args):
    if int(refGeneModelNo) == 3 and str(modelAnnoType).lower() == "gtf":
        return AGPV3_GTF
    elif int(refGeneModelNo) == 3 and str(modelAnnoType).lower() == "gff3":
        return AGPV3_GFF3
    elif int(refGeneModelNo) == 4 and str(modelAnnoType).lower() == "gtf":
        return AGPV4_GTF
    else:
        return AGPV4_GFF3

class GeneModel(object):
    def __init__(self):
        self.allGenes = {}
        self.allChr = {}
        self.geneModel = {}

    def updateChrom(self, start, end, chrId):
        chrId = chrId.lower()
        if chrId not in self.allChr:
            self.allChr[chrId] = Chromosome(start, end, chrId)
            self.geneModel.setdefault(chrId, {})
        else:
            record = self.allChr[chrId]
            record.minpos = min(record.minpos, start)
            record.maxpos = max(record.maxpos, end)
