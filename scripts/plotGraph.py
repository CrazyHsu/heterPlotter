#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
File name: plotGraph.py
Author: CrazyHsu @ crazyhsu9527@gmail.com 
Created on: 2018-10-16 17:32:40
Last modified: 2018-10-16 17:32:42
'''

from commonFunc import *
import pandas as pd
import numpy as np
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt

#================Variables=================
fontsize = 12
plotHeight, plotWidth = 11.0, 8.5

Y_LIMIT = 100.0
BORDER = 0.08
PLOT_AREA_SPECS = {'font.weight': 'normal', 'legend.borderaxespad': 0.01, 'legend.handlelength': 0.02, \
                     'legend.handletextpad': 0.01, 'legend.labelspacing': 0.008}
DISPLAY_HEIGHT = 0.98
AXIS_LEFT = 0.08
AXIS_WIDTH = 0.84
X_PAD_FRACTION = 0.01
Y_PAD_FRACTION = 0.0035

MIN_SEC_NUM = 5

#================Methods===================
def initPlots(width, height, **args):
    matplotlib.rcParams.update(PLOT_AREA_SPECS)
    matplotlib.rcParams['figure.figsize'] = width, height
    matplotlib.rcParams['font.size'] = 12
    matplotlib.rcParams['legend.fontsize'] = 3

def plotGeneDist(topLine, height, xBound, geneList):
    scaledLen = abs(xBound[1] - xBound[0])
    ax = plt.axes([AXIS_LEFT, topLine, AXIS_WIDTH, height])
    padding = X_PAD_FRACTION * scaledLen
    for gene in geneList:
        geneStart = gene.minpos
        geneEnd = gene.maxpos
        geneLen = geneEnd - geneStart
        arrowLen = geneLen/2.0
        if gene.strand == "+":
            ax.arrow(geneStart, 50, geneLen, 0.0, fc="black", ec="black", head_width=50,
                     head_length=arrowLen, length_includes_head=True, shape="full")
        else:
            ax.arrow(geneStart, 50, -geneLen, 0.0, fc="black", ec="black", head_width=50,
                     head_length=arrowLen, length_includes_head=True, shape="full")
    ax.set_xlim(xBound[0] - padding, xBound[1] + padding)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.xaxis.set_ticks_position("none")
    ax.set_ylim(0, 100)
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.yaxis.set_ticks_position("none")

# def plotGeneDist(topLine, height, xBound, geneList):
#     scaledLen = abs(xBound[1] - xBound[0])
#     ax = plt.axes([AXIS_LEFT, topLine, AXIS_WIDTH, height])
#     padding = X_PAD_FRACTION * scaledLen
#     for gene in geneList:
#         geneStart = gene.minpos
#         geneEnd = gene.maxpos
#         xStart = xBound[0]
#         xEnd = xBound[1]
#         scaledGeneStart = (geneStart - xStart) * float(scaledLen) / (xEnd - xStart)
#         scaledGeneEnd = (geneEnd - xStart) * float(scaledLen) / (xEnd - xStart)
#         scaledGeneLen = scaledGeneEnd - scaledGeneStart
#         arrowLen = scaledGeneLen/2.0
#         if gene.strand == "+":
#             ax.arrow(scaledGeneStart, 50, scaledGeneLen, 0.0, fc="black", ec="black", head_width=50,
#                      head_length=arrowLen, length_includes_head=True, shape="full")
#         else:
#             ax.arrow(scaledGeneEnd, 50, -scaledGeneLen, 0.0, fc="black", ec="black", head_width=50,
#                      head_length=arrowLen, length_includes_head=True, shape="full")
#     ax.set_xlim(0 - padding, scaledLen + padding)
#     ax.set_xticklabels([])
#     ax.set_xticks([])
#     ax.xaxis.set_ticks_position("none")
#     ax.set_ylim(0, 100)
#     ax.set_yticklabels([])
#     ax.set_yticks([])
#     ax.yaxis.set_ticks_position("none")

#================Classes===================
class PlotGraph(object):
    def __init__(self, hpFile, gwFile, geneList, xBound, topLine, height, geneTrackHeight, firstFlag, endFlag, **args):
        self.hpFile = hpFile
        self.gwFile = gwFile
        self.geneList = geneList
        self.xBound = xBound
        self.totalHeight = height
        self.topLine = topLine
        self.geneTrackHeight = geneTrackHeight
        self.firstFlag = firstFlag
        self.endFlag = endFlag
        self.heterHeight = self.totalHeight*4/10.0
        self.logHeight = self.totalHeight*6/10.0
        self.pmax = getAttribute("pmax", 15, **args)
        self.pmin = getAttribute("pmin", 5, **args)

    def plot(self):
        self.plotLogDist()
        self.plotHeterDist()

    def plotHeterDist(self):
        fileBase = os.path.basename(self.hpFile)
        fileWithoutExt = os.path.splitext(fileBase)[0]
        hp = pd.read_table(self.hpFile, index_col=0)
        pect = hp.iloc[:, -3:]
        pectValues = pect.values
        ax = plt.axes([AXIS_LEFT, self.topLine, AXIS_WIDTH, self.heterHeight])
        padding = X_PAD_FRACTION * len(pectValues)
        xvalues = setXticks(int(round(0 - padding)), int(round(len(pectValues) + padding)))
        legendLabels = pect.columns.values
        legendPatch = ax.stackplot(range(len(pectValues)), pectValues.T, linewidth=0)
        ax.set_xlim(0-padding, len(pectValues)-1+padding)
        ax.set_yticklabels([0, 1])
        ax.set_yticks([0, 1])
        ax.set_ylim(0, 1)
        if self.endFlag:
            start, end = self.xBound
            xvaluesReal = getXticksReal(start, end, xvalues)
            xlabels = ["%d" % x for x in xvaluesReal]
            plt.figlegend(legendPatch, legendLabels, loc=8, ncol=len(legendLabels), fontsize="medium",
                          bbox_to_anchor=(0.5, 0.01), handlelength=2, handletextpad=0.5,
                          fancybox=True, columnspacing=0.5)
            ax.set_title(fileWithoutExt, y=self.totalHeight / float(self.heterHeight))
            ax.set_xticklabels(xlabels)
            ax.set_xticks(xvalues)
        elif self.firstFlag:
            ax.set_title(fileWithoutExt, y=(self.totalHeight + self.geneTrackHeight) / float(self.heterHeight))
            ax.set_xticklabels([])
            ax.set_xticks([])
        else:
            ax.set_title(fileWithoutExt, y=self.totalHeight / float(self.heterHeight))
            ax.set_xticklabels([])
            ax.set_xticks([])

    def plotLogDist(self):
        gw = pd.read_table(self.gwFile, index_col=1)
        pos = gw.iloc[:, 1].values
        minpos = pos.min()
        maxpos = pos.max()
        p = gw.iloc[:, -1].values
        ppp = -np.log10(p)
        padding = X_PAD_FRACTION * (maxpos - minpos)
        ax = plt.axes([AXIS_LEFT, self.topLine + self.heterHeight, AXIS_WIDTH, self.logHeight])
        scatterSize = getScatterSize(len(ppp))
        ax.scatter(pos, ppp, s=scatterSize)
        ax.set_xlim(minpos - padding, maxpos + padding)
        ax.set_xticklabels([])
        ax.set_xticks([])
        ax.set_yticklabels([self.pmin + 2, self.pmax])
        ax.set_yticks([self.pmin + 2, self.pmax])
        ax.set_ylim(self.pmin, self.pmax)