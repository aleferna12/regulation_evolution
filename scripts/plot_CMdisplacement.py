#!/usr/bin/python2.7

'''
Plot the displacement of the total cell mass over a fixed period for one or more simulations.
simulations can be grouped by parameter set by numbering a subset of the files:
./plot_CMdisplacement.py 5000 1 file1 file2 file2 2 file3 file4 3 file5 file6 ....
'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np

#manager = plt.get_current_fig_manager()
#print manager
#sys.exit(1)
#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
fig, ax = plt.subplots()
if len(sys.argv) <3:
  print "This is the program 'plot_CMdisplacement.py'"
  print "Usage: ./plot_CMdisplacement.py <figure name> <timeperiod> <(short) simset <filename(s)>> "
  sys.exit(1)
else:
  figname=sys.argv[1]
  period=int(sys.argv[2])
  
labels=[]
displace=[]
offset=[]
colourset=[]
totalcounter=0
labelset=[]
for filename in sys.argv[3:]:

  if len(filename)<4:
    totalcounter+=1
    labelset.append(filename)
  else:
    #now we calculate the center of mass of all cells, and its displacement over x timesteps
    colourset.append(colours[totalcounter])

    offset.append(totalcounter+random.random()*0.6-0.3)
    startxpos=0.
    startypos=0.
    endxpos=0.
    endypos=0.
    
    count=0
    count2=0
    ##read data from file: store raw cell positions
    with open(filename,"r") as fin:

      #read first line first to get the first time point (could probably more clever, but hey)
      for line in fin:
        line=line.split()
        timepoint=int(line[0])
        startxpos+=float(line[3])
        startypos+=float(line[4])
        count+=1
        break

      #read rest of file, to find positions at beginning and end of period
      for line in fin:
        line=line.split()
        if int(line[0])==timepoint:
          startxpos+=float(line[3])
          startypos+=float(line[4])
          count+=1
        elif int(line[0])==timepoint+period:
          endxpos+=float(line[3])
          endypos+=float(line[4])
          count2+=1

    #calculate average positions at beginning and end of period
    startxpos/=float(count)
    startypos/=float(count)
    endxpos/=float(count2)
    endypos/=float(count2)

    xx=endxpos-startxpos
    yy=endypos-startypos
    
    displace.append(math.sqrt(xx*xx+yy*yy))
    
 
##start plotting##

ax.scatter(offset,displace, c=colourset,alpha=0.5)
ax.set_xlim([0,totalcounter+1])
#custom xlabels
my_labels=ax.get_xticks().tolist()
counter=0
newlab=[]
for el in my_labels:
  if el.is_integer() and counter<len(labelset) and el>0.001:
    newlab.append(labelset[counter])
    counter+=1
  else:
    newlab.append('')

ax.set_xticklabels(newlab)
#print my_labels
#print labelset
#print newlab

plt.xlabel('simulation')
ylab='displacement in time period '+str(period)
plt.ylabel (ylab)
fig.savefig(figname, bbox_inches='tight')
#plt.show()
