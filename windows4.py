#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 15:32:19 2018

@author: blaufer
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 13:37:51 2018

@author: lab

masking - https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.ma.masked_where.html#numpy.ma.masked_where

"""
import os
import csv
#import Set
import gzip
import itertools
import time
#import numpy.ma as ma
import multiprocessing as mp

startTime = time.time()

#wdir = '/Volumes/Shared/SOM/medmicro/labs/LaSalle Lab/Hyeyeon/Windows/hg38samples/'
wdir = os.getcwd() 
#wdir = '/Volumes/lasalle lab/Hyeyeon/Windows/hg38samples/'




#sampleSize = 4

# RAM per core 2-4 GB
sample = {'1':'JLTT0003_USPD16083426_HGJ52CCXY_L6_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.merged_CpG_evidence.cov.gz',
          #'2':'JLTT0001_USPD16083424_HGJ52CCXY_L6_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.merged_CpG_evidence.cov.gz',
          #gzip -cd JLTT0002_USPD16083425_HGJ52CCXY_L6_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.merged_CpG_evidence.cov.gz | head -n 966470 | gzip > J2_chr_90K.gz
          '2':'J2_chr_90K_USPD16083425_HGJ52CCXY_L6_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.merged_CpG_evidence.cov.gz',
          '3':'JLTT0002_USPD16083425_HGJ52CCXY_L6_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.merged_CpG_evidence.cov.gz',
          '4':'JLTT0004_USPD16083427_HGJ52CCXY_L6_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.merged_CpG_evidence.cov.gz'
         }
#       
#print sample['1']

cpiFile = wdir + 'hg38_genome_CGI.bed'
genome = 'hg38'
inFile = wdir + 'J2_chr_90K_USPD16083425_HGJ52CCXY_L6_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.merged_CpG_evidence.cov.gz'
     
# validChrm = chromosomes of interest
    # [ 1..22 or 1..20, X, Y, M ] 
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
    
# chromosome sizes used to partition windows
sizes = {}
with open(wdir + "hg38.chrom.sizes") as sizeFile:
    for line in sizeFile:

        if line.split('\t')[0].split('chr')[1] in validChrm:
            print line.split('\t')[0].split('chr')[1]
            (key, val) = line.split()
            partitions = int(val) / 20000
            sizes[key] = [int(val), partitions]
 
chrms = {}
chrmNumList = [] #unique valid chromosomes

with gzip.open(inFile, "rt") as bisInFile:
    reader = csv.reader(bisInFile, delimiter = '\t')
    
    for row in reader:
        chrmNum = row[0].split("chr")[1]   
        if chrmNum in validChrm:   
            validChrmNum = chrmNum   
            # add each row of a valid chromosome in the Bismark file to chrms[key = chrmNum]
            chrms.setdefault(validChrmNum,[]).append(row)
            if validChrmNum not in chrmNumList: #1st occurrence of chrm
                chrmNumList.append(validChrmNum)
     


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
            #cpis_Pos.setdefault(cpiLines,[]).append(range(int(cpiRow[1]), (int(cpiRow[2]) + 1)))
            #cpis_Pos.append(range(int(cpiRow[1]), (int(cpiRow[2]) + 1)))
    
    #print cpis       
    #print "cpis dictionary"
    #print cpis_Pos
#windows = {}
#windowCount = 0
#cpiL = 1
#for key in chrms: 
    #if windowCount < 2:      
    #for bRow in chrms[key]: #row in a specific chrm
        #if bRow[1] not in cpis_Pos.itervalues(): #row
            #sumReads = int(bRow[4]) + int(bRow[5])
            #if sumReads > 0:
                #windowCount = windowCount + 1
                #print windowCount
                #windows.setdefault(key,[]).append(float(bRow[3])/100)              
#print windows
#print cpis_Pos.itervalues()
#print cpis_Pos["15"]
#numList = sum(1 for x in cpis_Pos["15"] if isinstance(x, list))
#print numList
#print chrms["15"]
#bRowPos = []
rowNOcpi = []
rowYEScpi = []
print("15") 


flatCpiPos = list(itertools.chain.from_iterable(cpis_Pos["15"]))

setCpi = set(flatCpiPos)

print(flatCpiPos)
print(type(flatCpiPos[0]))


    
chromNo = "15"

######### test windows ##########
import gzip
import csv

inFile = wdir + 'J2_chr_90K_USPD16083425_HGJ52CCXY_L6_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.merged_CpG_evidence.cov.gz'

chrms = {}
chrmNumList = []
with gzip.open(inFile, "rt") as bisInFile:
    reader = csv.reader(bisInFile, delimiter = '\t')
    for row in reader:
        chrmNum = row[0].split("chr")[1]   
        if chrmNum in validChrm:   
            validChrmNum = chrmNum   
            chrms.setdefault(validChrmNum,[]).append(row)
            if validChrmNum not in chrmNumList: 
                chrmNumList.append(validChrmNum)

#wdir = '/Volumes/lasalle lab/Hyeyeon/Windows/hg38samples/'
# ********************* bounds ******************
import math
bounds = {}
for c in chrmNumList:
    firstChrmRow = chrms[c][0]
    lastChrmRow = chrms[c][-1]
    
    firstChrm = int(firstChrmRow[1])
    lastChrm = int(lastChrmRow[2])
    lowBound = int( math.floor(firstChrm / 20000) * 20000 )
    upBound = int( math.ceil(lastChrm / 20000) * 20000 )
    bounds[c] = [ lowBound, upBound ]

# ********************* write Header ************    
numSamples = 1  
with open(wdir + "chrmSizesTest.txt","wb") as outFile:
    writer = csv.writer(outFile, delimiter = '\t', lineterminator = '\n')
    headerData = [None] * (3 + numSamples)
    headerData[0] = 'chr'
    headerData[1] = 'start'
    headerData[2] = 'end'
    headerData[3] = 'JLTT0002'
    writer.writerow(headerData)
        
# ******************* write 20KB windows **************
with open(wdir + "chrmSizesTest.txt","ab") as outFile:
    writer = csv.writer(outFile, delimiter = '\t', lineterminator = '\n')
    for r in bounds.keys():  # r is the chromosome 
        
        for x in range(bounds[r][0], bounds[r][1]):

            outData = [None] * (3 + numSamples)
            outData[0] = r
            outData[1] = x
            outData[2] = 20000 + outData[1]
            writer.writerow(outData) 
            
# ************************ end of test windows ******************************
            
def checkCpi(bRow): 
    
    if int(bRow[1]) not in setCpi: 
        #print(bRow[1])
        rowNOcpi.append(bRow)
        print(rowNOcpi)
    elif int(bRow[1]) in setCpi:
        #print(bRow[1])
        rowYEScpi.append(bRow)
    
# write header for each file


    
cores = 12
pool = mp.Pool(processes=cores)
    
#pool.map(writeHeader, 'JLTT0002')
pool.map(checkCpi, chrms[chromNo])
pool.terminate()


print(rowNOcpi)
print(len(rowNOcpi))

print(rowYEScpi)
print(len(rowYEScpi))    

endTime = time.time()
print("Time = " + str(endTime-startTime))


    


                    
                    
           
