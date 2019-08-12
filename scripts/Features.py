#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
File name: Features.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2018-09-03 22:13:13
Last modified: 2018-09-03 22:13:14
'''
from commonFunc import getAttribute

class Chromosome(object):
    def __init__(self, start, end, chrId):
        self.minpos = start
        self.maxpos = end
        self.chrId = chrId

    def __len__(self):
        return abs(self.maxpos - self.minpos) + 1

    def __str__(self):
        return "{}:{}-{}".format(self.chrId, self.minpos, self.maxpos)

    def contains(self, pos):
        return self.minpos <= pos <= self.maxpos

class CDS(object):
    def __init__(self, start, end, chrId, strand, phase, attr={}):
        self.cdsAttr = attr
        self.chrId = chrId
        self.end = max(start, end)
        self.feature = "cds"
        self.parents = []
        self.phase = phase
        self.start = min(start, end)
        self.strand = strand

    def __str__(self):
        return "{} ({}):{}-{}({})".format(self.feature, self.chrId, self.start, self.end, self.strand)

    def __len__(self):
        return abs(self.end - self.start) + 1

class Exon(object):
    def __init__(self):
        self.childExons = {}
        self.child = self
        self.chrId = None
        self.end = 0
        self.exonAttr = {}
        self.feature = "exon"
        self.parentExons = {}
        self.parent = self
        self.parents = []
        self.start = 0
        self.strand = None

    def updateParents(self, trans):
        if trans not in self.parents:
            self.parents.append(trans)

    def __len__(self):
        return abs(self.start - self.end) + 1

    def __str__(self):
        return "{} ({}):{}-{}({})".format(self.feature, self.chrId, self.start, self.end, self.strand)

class ExonInfo(object):
    def __init__(self, exonList, **args):
        self.exons = self.getExonInfo(exonList, **args)
        self.minpos = min([n.start for n in self.exons])
        self.maxpos = max([n.end for n in self.exons])
        self.transId = getAttribute("transId", None, **args)
        self.strand = getAttribute("strand", ".", **args)
        self.splicePos = getAttribute("splicePos", None, **args)

    def getExonInfo(self, exonList, **args):
        myList = [Exon() for i in exonList]
        transId = getAttribute("transId", None, **args)
        if len(myList) != 1:
            for n in range(len(myList)):
                if n == 0:
                    myList[n].parent = myList[n]
                    myList[n].child = myList[n+1]
                    myList[n].parentExons[transId] = myList[n]
                    myList[n].childExons[transId] = myList[n+1]
                elif n == len(myList)-1:
                    myList[n].parentExons[transId] = myList[n-1]
                    myList[n].childExons[transId] = myList[n]
                    myList[n].parent = myList[n-1]
                    myList[n].child = myList[n]
                else:
                    myList[n].parentExons[transId] = myList[n-1]
                    myList[n].childExons[transId] = myList[n+1]
                    myList[n].parent = myList[n-1]
                    myList[n].child = myList[n+1]
                myList[n].start = exonList[n][0]
                myList[n].end = exonList[n][1]
        else:
            myList[0].parentExons[transId] = myList[0]
            myList[0].childExons[transId] = myList[0]
            myList[0].parent = myList[0]
            myList[0].child = myList[0]
            myList[0].start = exonList[0][0]
            myList[0].end = exonList[0][1]
        return myList

    @staticmethod
    def sortExon(exons, **args):
        transId = getAttribute("transId", None, **args)
        if len(exons) != 1:
            for n in range(len(exons)):
                if n == 0:
                    exons[n].parentExons[transId] = exons[n]
                    exons[n].childExons[transId] = exons[n+1]
                    exons[n].parent = exons[n]
                    exons[n].child = exons[n+1]
                elif n == len(exons) - 1:
                    exons[n].parentExons[transId] = exons[n-1]
                    exons[n].childExons[transId] = exons[n]
                    exons[n].parent = exons[n-1]
                    exons[n].child = exons[n]
                else:
                    exons[n].parentExons[transId] = exons[n-1]
                    exons[n].childExons[transId] = exons[n+1]
                    exons[n].parent = exons[n-1]
                    exons[n].child = exons[n+1]
        else:
            exons[0].parentExons[transId] = exons[0]
            exons[0].childExons[transId] = exons[0]
            exons[0].parent = exons[0]
            exons[0].child = exons[0]
        return exons

class Gene(object):
    def __init__(self, geneId, start, end, chrId, strand, geneName=None, attrs={}):
        self.geneId = geneId
        self.chrId = chrId
        self.cdsDict = {}
        self.minpos = min(start, end)
        self.maxpos = max(start, end)
        self.strand = strand
        self.attrs = attrs
        self.geneName = geneName
        self.trans = {}
        self.exons = []
        self.exonDict = {}
        self.cds = []
        self.utrs = []
        self.utrDict = {}
        self.startCodons = {}
        self.stopCodons = {}

    def updateTrans(self, transObj):
        transObj.parent = self
        return self.trans.setdefault(transObj.transId, transObj)

    def __str__(self):
        return "{} ({}):{}-{} (len={}, strand={})".format(self.geneId, self.chrId, self.minpos, self.maxpos, len(self), self.strand)

    def __len__(self):
        return abs(self.maxpos - self.minpos) + 1

class GeneInfo(object):
    def __init__(self, lineObj):
        self.chrId = lineObj.chrom
        self.geneName = lineObj.geneName
        # self.asType = lineObj.asType
        self.strand = lineObj.strand
        self.transList = lineObj.transList
        self.exonDict = {}
        self.minposList, self.maxposList = [], []

    def __len__(self):
        return abs(self.maxposList[-1] - self.minposList[0]) + 1

class Trans(object):
    def __init__(self, transId, start, end, chrId, strand, transAttr={}):
        self.transId = transId
        self.minpos = start
        self.maxpos = end
        self.chrId = chrId
        self.strand = strand
        self.transAttr = transAttr
        self.parent = None
        self.exons = []
        self.exonDict = {}
        self.utrs = []
        self.utrDict = {}
        self.cds = []
        self.cdsDict = {}
        self.startCodon = ()
        self.stopCodon = ()

    def updateExons(self, exonObj, **args):
        exonBound = (exonObj.start, exonObj.end)
        # print self.transId, self.exonDict
        transId = getAttribute("transId", None, **args)
        if exonBound not in self.exonDict:
            # preExon = self.exons[-1] if len(self.exons) != 0 else exonObj
            # preExon.child = exonObj
            # exonObj.parent = preExon
            self.exons.append(exonObj)
            self.exons = sorted(self.exons, key=lambda n: n.start)
            self.exonDict[exonBound] = exonObj
            if self not in exonObj.parents:
                exonObj.parents.append(self)
                self.minpos = min(self.minpos, exonObj.start)
                self.maxpos = max(self.maxpos, exonObj.end)
            # self.exons = ExonInfo.sortExon(self.exons)
            self.exons = ExonInfo.sortExon(self.exons, transId=transId)

    def updateCDS(self, cdsObj):
        cdsBound = (cdsObj.start, cdsObj.end)
        if cdsBound not in self.cdsDict:
            self.cdsDict[cdsBound] = cdsObj

        if self not in cdsObj.parents:
            cdsObj.parents.append(self)
        self.cds.append(cdsObj)
        self.cds = sorted(self.cds, key=lambda n: n.start)

    def updateUTR(self, utrObj):
        utrBound = (utrObj.start, utrObj.end)
        if utrBound not in self.utrDict:
            self.utrDict[utrBound] = utrObj

        if self not in utrObj.parents:
            utrObj.parents.append(self)
        self.utrs.append(utrObj)
        self.utrs = sorted(self.utrs, key=lambda n: n.start)

    def __len__(self):
        return abs(self.maxpos - self.minpos) + 1

    def __str__(self):
        return "{} ({}):{}-{}({}) len={}".format(self.transId, self.chrId, self.minpos, self.maxpos, self.strand, len(self))

class UTR(object):
    def __init__(self, start, end, chrId, strand, feature, attr={}):
        self.start = min(start, end)
        self.end = max(start, end)
        self.chrId = chrId
        self.strand = strand
        self.feature = feature
        self.utrAttr = attr
        self.parents = []

    def __len__(self):
        return abs(self.end - self.start) + 1

    def __str__(self):
        return "{}:{}-{} {}({})".format(self.chrId, self.start, self.end, self.feature, self.strand)