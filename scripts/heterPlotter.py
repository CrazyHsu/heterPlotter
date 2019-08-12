#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
File name: heterPlotter.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2018-10-15 13:51:01
Last modified: 2018-10-15 13:51:02
'''

import argparse, datetime
from GTFLoader import *
from plotGraph import *

#===============parse arguments=============
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gtf", default="gene.gtf", type=str,
                    help="This is the gtf file the species you specified.")
parser.add_argument("-i", "--interval", default=None, type=str,
                    help="This is the interval you specified.")
parser.add_argument("-a", "--anno", default="annotation.txt", type=str,
                    help="This is the annotation file you should provide.")
parser.add_argument("-idx", "--index", default=1000, type=int,
                    help="This is the index offset you want to build.")
parser.add_argument("-hp", "--hapmap", default="example.hapmap", type=str,
                    help="This is the config file which contains the absolute path of hapmap files.")
parser.add_argument("-gw", "--gwas", default="example.gwas", type=str,
                    help="This is the config file which contains the absolute path of gwas files.")
parser.add_argument("-pmax", "--pvalueMax", default=15, type=int,
                    help="This is the max pvalue you specified to draw the plot. Default is 15.")
parser.add_argument("-pmin", "--pvalueMin", default=5, type=int,
                    help="This is the min pvalue you specified to draw the plot. Default is 5.")
parser.add_argument("-o", "--out", default="figure.png", type=str,
                    help="This is the output figure you specified.")

#=================Main Methods==============
def main():
    startTime = datetime.datetime.now()
    args = parser.parse_args()
    assert validateFile(args.gtf) and validateFile(args.hapmap) and validateFile(args.gwas)
    hapmapList = getFileList(args.hapmap)
    gwasList = getFileList(args.gwas)
    pmax, pmin = args.pvalueMax, args.pvalueMin
    idxNum = args.index
    assert int(pmax) and int(pmin)
    gtfStartTime = datetime.datetime.now()
    gtfModel = GTFLoader(args.gtf)
    gtfEndTime = datetime.datetime.now()
    print "GTF loader runs %s" % (gtfEndTime - gtfStartTime).total_seconds() + " seconds."
    interval = args.interval
    checkInterval(interval, gtfModel)
    chrId, start, end = getPos(interval)
    curDir = os.getcwd()
    annoFile = args.anno

    hpTmpPectList, tmpGwasList = [], []
    maxPosList, minPosList = [], []
    genoTypeDict = {}
    tmpStartTime = datetime.datetime.now()
    for i in range(len(hapmapList)):
        hapmap = hapmapList[i]
        assert validateFile(hapmap)
        hapmapBase = os.path.basename(hapmap)
        hpIdxFile = os.path.join(curDir, hapmapBase + ".idx")
        hpFileTmp = os.path.join(curDir, hapmapBase + ".tmp")
        hpFileTmpPect = os.path.join(curDir, hapmapBase + ".pect")
        hpLine2Pos = buildIndex(hapmap, hpIdxFile, fileType=1, lineOffset=idxNum)
        if i == 0:
            maxPosTmp, minPosTmp, genoTypeDict = createTmpFile(interval, hpIdxFile, hapmap, hpFileTmp, fileType=1,
                                                               genoTypeDict={}, line2Pos=hpLine2Pos)
        else:
            maxPosTmp, minPosTmp, tmpDict = createTmpFile(interval, hpIdxFile, hapmap, hpFileTmp, fileType=1,
                                                          genoTypeDict=genoTypeDict, line2Pos=hpLine2Pos)
        hpTmpPectList.append(hpFileTmpPect)
        maxPosList.append(maxPosTmp)
        minPosList.append(minPosTmp)
    for gwas in gwasList:
        gwasBase = os.path.basename(gwas)
        assert validateFile(gwas)
        gwIdxFile = os.path.join(curDir, gwasBase + ".idx")
        gwFileTmp = os.path.join(curDir, gwasBase + ".tmp")
        gwLine2Pos = buildIndex(gwas, gwIdxFile, fileType=2, lineOffset=idxNum)
        maxPosTmp, minPosTmp = createTmpFile(interval, gwIdxFile, gwas, gwFileTmp, fileType=2, lineOffset=idxNum,
                                             line2Pos=gwLine2Pos)
        tmpGwasList.append(gwFileTmp)
        maxPosList.append(maxPosTmp)
        minPosList.append(minPosTmp)

    initPlots(plotWidth, plotHeight)
    height = DISPLAY_HEIGHT * (1.0/(len(hpTmpPectList)+2)) if len(hpTmpPectList) > MIN_SEC_NUM \
        else DISPLAY_HEIGHT * (1.0/(MIN_SEC_NUM+2))
    tmpEndTime = datetime.datetime.now()
    print "Kinds of temporal files generating expends %s" % (tmpEndTime - tmpStartTime).total_seconds() + " seconds."

    headPadding = height/(MIN_SEC_NUM-1)
    maxPos = max(maxPosList)
    minPos = min(minPosList)
    newInterval = chrId + ":" + str(minPos) + "-" + str(maxPos)
    geneList = getTargetGenesByPos(newInterval, gtfModel)
    geneList = sortGeneList(geneList)
    xBound = (minPos, maxPos)

    plt.figure()
    geneTrackHeight = float(height)/10
    topLine = DISPLAY_HEIGHT - headPadding - geneTrackHeight
    plotGeneDist(topLine, geneTrackHeight, xBound, geneList)
    topLine = topLine - height
    for i in range(len(hpTmpPectList)):
        hpFile = hpTmpPectList[i]
        gwFile = tmpGwasList[i]
        firstFlag = 1 if i == 0 else 0
        endFlag = 1 if i == len(hpTmpPectList) - 1 else 0
        pg = PlotGraph(hpFile, gwFile, geneList, xBound, topLine, height, geneTrackHeight, firstFlag, endFlag,
                       pmax=pmax, pmin=pmin)
        pg.plot()
        topLine = topLine - height - headPadding
    plt.savefig(args.out)

    if validateFile(annoFile):
        annoDict = getAnnoDict(annoFile)
        geneOut = open("targetGene.anno", "wb")
        for gene in geneList:
            if gene.geneName in annoDict:
                print >>geneOut, annoDict[gene.geneName]
        geneOut.close()
    endTime = datetime.datetime.now()
    print "This program runs %s" % (endTime-startTime).total_seconds() + " seconds."

if __name__ == "__main__":
    main()
