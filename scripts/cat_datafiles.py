#!/usr/bin/python2.7

'''
Script to merge multiple datafiles for concurrent analysis (general stats like MSD, for which sigma is not used)
Best to use two files from runs with the same parameters, otherwise it makes little sense.
Probably also better if they start at the same time.

'''

import sys,math,os,subprocess,random
#from PIL import Image
import numpy as np

outputfile=open(sys.argv[1],"w")

lastpoint=99999999
timepoint=0
duration=0
#find the timepoint that all files still have
for filename in sys.argv[2:]:

  with open(filename, "r") as fin:
    thisfile=list(fin)
    firstline=thisfile[0].split()
    lastline=thisfile[-1].split()

    if int(lastline[0])<lastpoint:
      lastpoint=int(lastline[0])
    if int(firstline[0])>timepoint:
      timepoint=int(firstline[0])
    if(not duration):
      secline=thisfile[1].split()
      duration=int(secline[0])-int(firstline[0])

print "first and last time point are ",timepoint, " ", lastpoint

#concatenate the files by timepoint
nexttime=timepoint
cellcount=1
prevcell=0
while timepoint<=lastpoint:
  #print timepoint
  for filename in sys.argv[2:]:
    #print filename
    with open(filename, "r") as fin:
      for line in fin:
        larr=line.split()
        if int(larr[0])==timepoint:
          if(int(larr[1])-prevcell>1):
            cellcount+=1
          prevcell=int(larr[1])
          larr[1]=str(cellcount)
          #print larr
          outputfile.write(' '.join(larr))
          outputfile.write("\n")
          cellcount+=1
        elif int(larr[0])>timepoint:
          nexttime=int(larr[0])
          #print "found nexttime in file ",filename,"it's ",nexttime
          break	
  
  timepoint+=duration
  cellcount=1
  prevcell=0
  #print timepoint

outputfile.close()
