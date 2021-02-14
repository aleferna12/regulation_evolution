#!/usr/bin/python2.7

'''
Script to merge multiple datafiles for concurrent analysis (general stats like MSD, for which sigma is not used)
Best to use two files from runs with the same parameters, otherwise it makes little sense.
Probably also better if they start at the same time.

'''

import sys,math,os,subprocess,random
#from PIL import Image
import numpy as np

inputfile=sys.argv[1]
#outputfile=open(sys.argv[2],"w")
season_start=int(sys.argv[2])
season_end=int(sys.argv[3])

with open(inputfile, "r") as fin:
  for line in fin:
    larr=line.split()
    if int(larr[0])<season_start:
      continue
    elif int(larr[0])>=season_start and int(larr[0])<season_end:
      print line,
    else:
      break
     
