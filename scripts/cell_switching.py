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
outputformat=int(sys.argv[4])

if (outputformat!=1 and outputformat!=2):
  print "wrong outputformat"
  sys.exit(1)


curr_cell=0
cellstor=[[] for x in range(5000)]
celltime=[0 for x in range(5000)]
cellswitch=[[] for x in range(5000)]
cellstart=[-1 for x in range(5000)]
with open(inputfile, "r") as fin:
  for line in fin:
    larr=line.split()
    if int(larr[0])<season_start:
      continue
    elif int(larr[0])>=season_start and int(larr[0])<season_end:
      #print larr[0]
      sigma=int(larr[1])
      tau=int(larr[2])
      cellstor[sigma].append(tau)
      celltime[sigma]+=1
      if cellstart[sigma]==-1:
        cellstart[sigma]=int(larr[0])-season_start
      if len(cellstor[sigma])>1:
	if cellstor[sigma][-1]!=cellstor[sigma][-2]:
          cellswitch[sigma].append(celltime[sigma]*(cellstor[sigma][-1]-cellstor[sigma][-2]))
          celltime[sigma]=0
        
    else:
      for sigma,taus in enumerate(cellstor):
        if len(taus):
          cellswitch[sigma].append(celltime[sigma]*(-1*(-1)**taus[-1] )  )
      break
      
if(outputformat==1):
  for ind,cel in enumerate(cellswitch):
    if len(cel)>1:
      print ind, cellstart[ind], cel
elif(outputformat==2):
  for ind,cel in enumerate(cellswitch):
    if len(cel)>0:
      for el in cel:
        print el 
