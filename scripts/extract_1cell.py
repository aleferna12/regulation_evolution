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
outputfile=open(sys.argv[2],"w")
sigma=int(sys.argv[3])
#newid=sys.argv[4]

with open(inputfile, "r") as fin:
  for line in fin:
    larr=line.split()
    if int(larr[1])==sigma:
      larr[1]="1"
      outputfile.write(' '.join(larr)+"\n")


outputfile.close()
