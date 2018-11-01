#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 09:50:01 2018

@author: lab
"""

import os
import csv
#import Set
import gzip
import itertools
import time
#import numpy.ma as ma
import multiprocessing as mp
import math

wdir = os.getcwd() 
#path = '/Volumes/lasalle lab/Hyeyeon/Windows/'

cpiFile = wdir + '/' + 'hg38_genome_CGI.bed'
genome = 'hg38'
inFile = wdir + '/' + 'J2_415.txt'

# validChrm = chromosomes of interest
validChrm = []
if genome == "hg38" or genome == "hg19":
    for i in range(1, 23):  #1..22
        validChrm.append(str(i))
if genome == "rheMac8" or genome == "rn6":
    for i in range(1, 21): #1..20
        validChrm.append(str(i))
if genome == "mm9" or genome == "mm10":
    for i in range(1, 20): #1..19
        validChrm.append(str(i))        
validChrm.append("X") 
validChrm.append("Y")
validChrm.append("M")

chrms = {}
chrmNumList = [] #unique valid chromosomes
with open(inFile, "rt") as bisInFile:
    reader = csv.reader(bisInFile, delimiter = '\t')
    
    for row in reader:
        chrmNum = row[0].split("chr")[1]   
        if chrmNum in validChrm:   
            validChrmNum = chrmNum   
            # add each row of a valid chromosome in the Bismark file to chrms[key = chrmNum]
            chrms.setdefault(validChrmNum,[]).append(row)
            if validChrmNum not in chrmNumList: #1st occurrence of chrm
                chrmNumList.append(validChrmNum)
               
# get CpG island ranges to mask                
with open(cpiFile, "rt") as cpiInFile:
    cpis_Pos = {}
    cpiLines = 0
    cpis = {}
    cpiReader = csv.reader(cpiInFile, delimiter = '\t') 
    for cpiRow in cpiReader:
        
        if cpiRow[0].split("chr")[1] in chrms.keys():
            cpiLines = cpiLines + 1
            validCpiNum = cpiRow[0].split("chr")[1]
            cpis.setdefault(validCpiNum,[]).append(cpiRow)
            cpis_Pos.setdefault(validCpiNum,[]).append(range(int(cpiRow[1]), (int(cpiRow[2]) + 1)))
            
rowNOcpi = []
rowYEScpi = []
print("15") 

flatCpiPos = list(itertools.chain.from_iterable(cpis_Pos["15"]))

setCpi = set(flatCpiPos)

print(flatCpiPos)
print(type(flatCpiPos[0]))

# ********************* bounds ******************

bounds = {}
for c in chrmNumList:
    firstChrmRow = chrms[c][0]
    lastChrmRow = chrms[c][-1]
    
    firstChrm = int(firstChrmRow[1])
    lastChrm = int(lastChrmRow[2])
    lowBound = int( math.floor(firstChrm / 20000.0) * 20000.0 )
    upBound = int( math.ceil(lastChrm / 20000.0) * 20000.0 )
    bounds[c] = [ lowBound, upBound ] 

# ********************* write Header ************    
numSamples = 1  
with open("windows-415-test.txt","wb") as outFile:
    writer = csv.writer(outFile, delimiter = '\t', lineterminator = '\n')
    headerData = [None] * (3 + numSamples)
    headerData[0] = 'chr'
    headerData[1] = 'start'
    headerData[2] = 'end'
    headerData[3] = 'JLTT0002'
    writer.writerow(headerData)
        
# ******************* write 20KB windows **************
outData415 = {} 
twentyWindow = {}
#outData4 = {}
with open("windows-415-test.txt","ab") as outFile:
    writer = csv.writer(outFile, delimiter = '\t', lineterminator = '\n')
    # dictionary for each chromosome at end?
    # dict for each bound range?
    
    # will contain bound range as keys?
    # add everything 
    increment = 0
    for r in bounds.keys(): # not bounds.keys, r is the chromosome 
        print bounds
        
        for x in xrange(bounds[r][0], bounds[r][1], 20000):
            outData = [None] * (3 + numSamples)
            outData[0] = r
            outData[1] = x
            outData[2] = x + 20000
            
            twentyWindow.setdefault(r, []).append([outData[1], outData[2]])
         
            outData415.setdefault(r +'-'+ str(outData[1]) +'-'+ str(outData[2]), [])
            
            increment = increment + 1
            
            writer.writerow(outData) 
         
for key in twentyWindow:
    for tRow in twentyWindow[key]:
        #print tRow
        #print str(tRow[1]) +'-'+ str(tRow[2])
        
        for chRow in chrms[key]:
            #print chRow
           
            if int(chRow[1]) in range(tRow[0], tRow[1]):
                keyName = key +'-'+ str(tRow[0]) +'-'+ str(tRow[1]) 
                outData415[keyName].append(chRow)
                #print chRow    
                #print tRow
            


         
        #if chrRow[1] in bounds, put it in that key of the outData415 dictionary
                     





