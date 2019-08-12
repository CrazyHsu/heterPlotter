#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
File name: GTF_Rec.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2018-09-03 21:42:31
Last modified: 2018-09-03 21:43:34
'''

GENE_ID_ATTR, TRANS_ID_ATTR, EXON_ID_ATTR, PROTEIN_ID_ATTR = \
    "gene_id", "transcript_id", "exon_id", "protein_id"
EXON_NO_ATTR, GENE_NAME_ATTR, TRANS_NAME_ATTR = \
    "exon_number", "gene_name", "transcript_name"
GENE_BIOTYPE, TRANS_BIOTYPE = "gene_biotype", "transcript_biotype"
GENE_SOURCE, TRANS_SOURCE = "gene_source", "transcript_source"

ID_ATTR = [GENE_ID_ATTR, TRANS_ID_ATTR]
NAME_ATTR = [GENE_NAME_ATTR, TRANS_NAME_ATTR]
GENE_OPT_ATTR = [GENE_NAME_ATTR, GENE_BIOTYPE, GENE_SOURCE]
TRANS_OPT_ATTR = [TRANS_NAME_ATTR, TRANS_BIOTYPE, TRANS_SOURCE]
NAME2ID = {GENE_NAME_ATTR: GENE_ID_ATTR, TRANS_NAME_ATTR: TRANS_ID_ATTR}

class GTF_Line(object):
    def __init__(self, line):
        if line.startswith("#"):
            raise ValueError()
        self.fields = line.strip().split("\t")
        if len(self.fields) != 9:
            raise ValueError("This is incorrect GTF format where the correct GTF format should contain 9 column")
        self.setFieldInfo()
        self.attrs = self.setAttrs()
        self.setAttrInfo()

    def setAttrs(self):
        attrList = self.fields[8].split(";")
        attrDict = {}
        for p in attrList:
            if not p: continue
            pair = p.split()
            key = pair[0].strip()
            value = pair[1].strip().replace('"', '')
            attrDict[key] = value.upper()

        absentList = []
        for attr in ID_ATTR:
            if attr not in attrDict:
                if self.feature.upper() == "GENE":
                    attrDict["transcript_id"] = None
                else:
                    absentList.append(attr)
        if len(absentList) > 0:
            raise ValueError("Please check if you have input the right GTF file with {} attribute".format(", ".join(absentList)))

        for attr in NAME_ATTR:
            if attr not in attrDict:
                attrDict[attr] = attrDict[NAME2ID[attr]]
        return attrDict

    def setFieldInfo(self):
        self.chrId = self.fields[0]
        self.source = self.fields[1]
        self.feature = self.fields[2]
        self.start = int(self.fields[3])
        self.end = int(self.fields[4])
        self.score = self.fields[5]
        self.strand = self.fields[6]
        self.phase = self.fields[7]

    def setAttrInfo(self):
        self.geneId = self.attrs[GENE_ID_ATTR]
        self.transId = self.attrs[TRANS_ID_ATTR]