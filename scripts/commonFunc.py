#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
File name: commonFunc.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2018-10-15 14:03:48
Last modified: 2018-10-15 14:03:52
'''

import os, re
import pandas as pd
from collections import Counter

#==================Variables================
HETERDICT={"AC": 1, "AG": 1, "AT": 1, "CA": 1, "CG": 1, "CT": 1,
           "TA": 1, "TC": 1, "TG": 1, "GA": 1, "GC": 1, "GT": 1}
LINEOFFSET = 1000
POINTSIZE = {2000: 2, 5000: 1.5, 8000: 1}

#==================Methods==================

def buildIndex(hapmapFile, out, fileType=1, lineOffset=LINEOFFSET):
    fileObj = open(hapmapFile, "rb")
    dirName = os.path.dirname(hapmapFile)
    outFile = os.path.join(dirName, out)
    fileOut = open(outFile, "wb")
    offsetList = []
    offset = 0
    chrFlag = True
    preChr = ""
    lineCount = 0
    chrPointList = []
    line2Pos = []
    for line in fileObj:
        offsetList.append(offset)
        offset += len(line)
        if fileType != 1:
            # tmpChr = line.split("\t")[0]
            infoList = line.strip("\n").split("\t")
            tmpChr = infoList[0]
            pos = infoList[2]
            line2Pos.append(pos)
            if chrFlag:
                preChr = tmpChr
                chrFlag = False
                lineCount += 1
                continue
            if tmpChr != preChr:
                preChr = tmpChr
                chrPointList.append(lineCount)
        else:
            infoList = line.strip("\n").split("\t")
            pos = infoList[3]
            line2Pos.append(pos)
        lineCount += 1
    fileObj.seek(0)

    offsetStart, numCount = 0, 0
    flag = True
    popFlag, firstLine = False, True
    if fileType == 1:
        while flag:
            try:
                if firstLine:
                    firstLine = False
                    offsetStart = 1
                    fileObj.seek(offsetList[offsetStart])
                    myLine = fileObj.readline().strip()
                    chrId = myLine.split("\t")[2]
                    pos = myLine.split("\t")[3]
                    print >> fileOut, "\t".join(
                        [str(i) for i in [numCount, offsetStart, pos, chrId, offsetList[offsetStart]]])
                    offsetStart = offsetStart + lineOffset
                    numCount += 1
                    continue
                fileObj.seek(offsetList[offsetStart])
                myLine = fileObj.readline().strip()
            except IndexError:
                fileObj.seek(offsetList[-1])
                myLine = fileObj.readline().strip()
                offsetStart = len(offsetList) - 1
                flag = False
            chrId = myLine.split("\t")[2]
            pos = myLine.split("\t")[3]
            print >>fileOut, "\t".join([str(i) for i in [numCount, offsetStart, pos, chrId, offsetList[offsetStart]]])
            offsetStart += lineOffset
            numCount += 1
    else:
        changeLine = chrPointList[0]
        chrPointList.pop(0)
        while flag:
            if popFlag and chrPointList:
                changeLine = chrPointList[0]
                chrPointList.pop(0)
                popFlag = False
            try:
                if offsetStart >= changeLine and chrPointList:
                    popFlag = True
                    fileObj.seek(offsetList[changeLine])
                    myLine = fileObj.readline().strip()
                    chrId = myLine.split("\t")[0]
                    pos = myLine.split("\t")[2]
                    print >> fileOut, "\t".join(
                        [str(i) for i in [numCount, changeLine, pos, chrId, offsetList[changeLine]]])
                    offsetStart = changeLine + lineOffset
                    numCount += 1
                    continue
                if firstLine:
                    firstLine = False
                    offsetStart = 1
                    # fileObj.seek(offsetList[offsetStart])
                    # myLine = fileObj.readline().strip()
                    # chrId = myLine.split("\t")[0]
                    # pos = myLine.split("\t")[2]
                    # print >> fileOut, "\t".join(
                    #     [str(i) for i in [numCount, offsetStart, pos, chrId, offsetList[offsetStart]]])
                    # offsetStart = offsetStart + lineOffset - 1
                    numCount += 1
                    continue
                fileObj.seek(offsetList[offsetStart])
                myLine = fileObj.readline().strip()
            except IndexError:
                fileObj.seek(offsetList[-1])
                myLine = fileObj.readline().strip()
                offsetStart = len(offsetList) - 1
                flag = False
            chrId = myLine.split("\t")[0]
            pos = myLine.split("\t")[2]
            print >>fileOut, "\t".join([str(i) for i in [numCount, offsetStart, pos, chrId, offsetList[offsetStart]]])
            offsetStart += lineOffset
            numCount += 1
    fileObj.close()
    fileOut.close()
    return line2Pos

def checkInterval(interval, gtfModel):
    interList = re.split("[:|-]", interval)
    chrId, start, end = interList[0], int(interList[1]), int(interList[1])
    allChr = gtfModel.model.allChr
    if chrId not in allChr:
        raise Exception("You should input the right chromosome id you specified!")
    else:
        if start < allChr[chrId].minpos or end > allChr[chrId].maxpos:
            raise Exception("You should input the right position you specified!")

def calMajorAndMinorAllele(x, genoTypeDict, femaleFlag=0):
    name = x.name
    counter = Counter(x).most_common(3)
    if femaleFlag:
        majorGeno, minorGeno, majorCount, minorCount = judgeMajorAndMinorParent(counter, len(counter))
        genoTypeDict[name] = majorGeno
    else:
        femaleMajorGeno = genoTypeDict[name]
        majorGeno, minorGeno, majorCount, minorCount = judgeMajorAndMinorChild(counter, len(counter), femaleMajorGeno)
    majorPect = round(majorCount / float(len(x)), 4)
    minorPect = round(minorCount / float(len(x)), 4)
    heterPect = round(1 - majorPect - minorPect, 4)
    return ",".join([majorGeno, minorGeno, str(majorPect), str(minorPect), str(heterPect)])

def convertHapmap2Num(hapmapFile, hapmap2Num, pectFile, genoTypeDict):
    data = pd.read_table(hapmapFile, index_col=0)
    pos = data.iloc[:, 2]
    headerdf = data.iloc[:, 0:10]
    genoType = data.iloc[:, 10:]
    probeId = data.index
    femaleFlag = 1 if len(genoTypeDict) == 0 else 0
    tmp = genoType.apply(lambda x: calMajorAndMinorAllele(x, genoTypeDict, femaleFlag=femaleFlag), axis=1)
    tmpdf = tmp.str.split(",", expand=True)
    tmpdf = tmpdf.rename({0: "majorGeno", 1: "minorGeno", 2: "majorPect", 3: "minorPect", 4: "heterPect"},
                         axis="columns")
    genoType = pd.concat([genoType, tmpdf], axis=1)
    grouped = genoType.groupby(["majorGeno", "minorGeno"])
    newdf = pd.DataFrame()
    for name, group in grouped:
        majorGeno = group["majorGeno"].unique()[0]
        minorGeno = group["minorGeno"].unique()[0]
        genoDict = {majorGeno: 0, minorGeno: 2}
        genoDictCopy = genoDict.copy()
        genoDictCopy.update(HETERDICT)
        newGroup = group.iloc[:, :-5].replace(genoDictCopy)
        newdf = newdf.append(newGroup)
    newdf = newdf.loc[probeId]
    finalGenodf = pd.concat([headerdf, newdf, tmpdf], axis=1)
    finalGenodf.to_csv(hapmap2Num, sep="\t")
    tmpdf.insert(loc=2, column=pos.name, value=pos)
    tmpdf.iloc[:, -4:].to_csv(pectFile, sep="\t")
    return genoTypeDict

def createTmpFile(interval, idxFile, originFile, fileTmp, **args):
    fileType = getAttribute("fileType", 1, **args)
    genoTypeDict = getAttribute("genoTypeDict", {}, **args)
    line2Pos = getAttribute("line2Pos", {}, **args)
    chrId, start, end = getPos(interval)
    tmpOut = open(fileTmp, "wb")
    curDir = os.getcwd()
    originFileBase = os.path.basename(originFile)
    numFile = os.path.join(curDir, originFileBase+".num")
    pectFile = os.path.join(curDir, originFileBase+".pect")
    preOffset, postOffset = 0, 0
    maxPos, minPos = 0, float("inf")
    with open(idxFile, "rb") as f1:
        chrFlag = 0
        for line in f1.readlines():
            sList = re.split("\t", line.strip())
            num, tmpLineIdx, pos, tmpChrId, tmpOffset = sList
            pos = int(pos)
            tmpOffset = int(tmpOffset)
            if tmpChrId == str(chrId):
                chrFlag = 1
                if int(pos) < start:
                    preOffset = tmpOffset
                    minPos = pos
                    continue
                if int(pos) < end:
                    postOffset = tmpOffset
                    maxPos = int(line2Pos[int(tmpLineIdx) - 1])
                    continue
                else:
                    maxPos = int(line2Pos[int(tmpLineIdx) - 1])
                    postOffset = tmpOffset
                    break
        if not chrFlag:
            raise Exception("The chromosome id you input doesn't contain in the %s" % originFile)
    with open(originFile, "rb") as f2:
        print >>tmpOut, f2.readline().strip()
        f2.seek(0)
        f2.seek(preOffset)
        print >>tmpOut, f2.read(postOffset-preOffset-1)
    tmpOut.close()
    if fileType == 1:
        if len(genoTypeDict) == 0:
            genoTypeDict = convertHapmap2Num(fileTmp, numFile, pectFile, genoTypeDict)
            return maxPos, minPos, genoTypeDict
        else:
            tmpDict = convertHapmap2Num(fileTmp, numFile, pectFile, genoTypeDict)
            return maxPos, minPos, tmpDict
    else:
        return maxPos, minPos

def getAnnoDict(annoFile):
    with open(annoFile) as f:
        annoDict = {}
        for i in f:
            infoList = re.split("\t", i.strip())
            if infoList[1] == "NA": continue
            annoDict[infoList[1]] = i.strip()
        return annoDict

def getAttribute(key, default, **args):
    return default if key not in args else args[key]

def getPos(interval):
    interList = re.split("[:|-]", interval)
    chrId, start, end = interList[0], interList[1], interList[2]
    return chrId, int(start), int(end)

def getScatterSize(pointNum):
    defaultSize = 0.5
    if pointNum < 2000:
        pointSize = POINTSIZE[2000]
    elif pointNum < 5000:
        pointSize = POINTSIZE[5000]
    elif pointNum < 8000:
        pointSize = POINTSIZE[8000]
    else:
        pointSize = defaultSize
    return pointSize

def getTargetGenesByPos(interval, gtfModel):
    chrId, start, end = getPos(interval)
    genesOnTargetChr = gtfModel.model.geneModel[chrId]
    specBound = (start, end)
    targetGenes = []
    for gene in genesOnTargetChr:
        geneBound = (genesOnTargetChr[gene].minpos, genesOnTargetChr[gene].maxpos)
        if judgeOverlap(specBound, geneBound):
            targetGenes.append(genesOnTargetChr[gene])
    return targetGenes

def getXticksReal(start, end, xvalues):
    xvaluesReal = []
    for i in range(len(xvalues)):
        xReal = (end - start) * i/(len(xvalues) - 1) + start
        xvaluesReal.append(xReal)
    return xvaluesReal

def judgeOverlap(tupleA, tupleB):
    return max(0, min(tupleA[1], tupleB[1]) - max(tupleA[0], tupleB [0]))

def judgeMajorAndMinor(counter, counterLen):
    returnList = []
    if counterLen == 3:
        if counter[0][0] in HETERDICT:
            returnList = [counter[1][0], counter[2][0], counter[1][1], counter[2][1]]
        elif counter[1][0] in HETERDICT:
            returnList = [counter[0][0], counter[2][0], counter[0][1], counter[2][1]]
        elif counter[2][0] in HETERDICT:
            returnList = [counter[0][0], counter[1][0], counter[0][1], counter[1][1]]
    elif counterLen == 2:
        if counter[0][0] in HETERDICT:
            minorAllele = ''
            for i in counter[0][0]:
                if i not in counter[1][0]:
                    minorAllele = i
            returnList = [counter[1][0], minorAllele*2, counter[1][1], 0]
        elif counter[1][0] in HETERDICT:
            minorAllele = ''
            for i in counter[1][0]:
                if i not in counter[0][0]:
                    minorAllele = i
            returnList = [counter[0][0], minorAllele * 2, counter[0][1], 0]
        else:
            returnList = [counter[0][0], counter[1][0], counter[0][1], counter[1][1]]
    elif counterLen == 1:
        returnList = [counter[0][0], "NN", counter[0][1], 0]
    return returnList

def judgeMajorAndMinorParent(counter, counterLen):
    return judgeMajorAndMinor(counter, counterLen)

def judgeMajorAndMinorChild(counter, counterLen, femaleMajorGeno):
    majorGeno, minorGeno, majorCount, minorCount = judgeMajorAndMinor(counter, counterLen)
    if majorGeno == femaleMajorGeno:
        return majorGeno, minorGeno, majorCount, minorCount
    elif minorGeno == femaleMajorGeno:
        return minorGeno, majorGeno, minorCount, majorCount

def setXticks(minpos, maxpos, minTicks=4, maxTicks=10) :
    """
    Overrides the matplotlib default to use absolute (rather than relative) tick
    mark labels when positions exceed ~10^7
    """
    deltas = [1000000, 500000, 100000, 50000, 10000, 5000, 2000, 1000, 500, 200, 100, 50, 20, 10]
    result = range(minpos,maxpos,10)
    for d in deltas :
        firstTick = d + d*int(int(minpos)/d)
        result = range(firstTick,int(maxpos)+1,d)
        if len(result) >= minTicks : break

    # Remove every other tick label
    if len(result) > maxTicks :
        revised = [result[i] for i in xrange(0,len(result),2)]
        result = revised
    return result

def sortGeneList(geneObjList):
    return sorted(geneObjList, key=lambda x: x.minpos)

def validateFile(myFile):
    if not os.path.exists(myFile):
        raise Exception("File '%s' not found! Please input again!" % myFile)

    if not os.path.isfile(myFile):
        raise Exception("File '%s' is not a file! Please input again!" % myFile)

    return True

def validateDir(myDir):
    if not os.path.exists(myDir):
        raise Exception("Dir '%s' not found! Please input again!" % myDir)

    if not os.path.isdir(myDir):
        raise Exception("Dir '%s' is not a directory! Please input again!" % myDir)

    return True

def getFileList(myFile):
    with open(myFile) as f:
        myList = []
        for line in f:
            if line.startswith("#"): continue
            assert validateFile(line.strip())
            myList.append(line.strip())
        return myList
