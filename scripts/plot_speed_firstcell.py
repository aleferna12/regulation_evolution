#!/usr/bin/python2.7

'''
Plot the instant speed fluctuations
Plot some selected cell tracks, centered at 0,0
https://www.itl.nist.gov/div898/handbook/eda/section3/eda35c.htm
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

colours=[(0.7,0.13,0.),"royalblue", "darkgoldenrod", (0.,0.5,0.2), "salmon", "lightskyblue","orchid","chocolate"]
filename=""
fig, (ax0, ax1) = plt.subplots(ncols=2,sharey=True,sharex=True)
fig.subplots_adjust(wspace=0)
#fig, ax0 = plt.subplots()

if len(sys.argv) <3:
  print "This is the program 'plot_speed.py'"
  print "Usage: ./plot_speed.py <figure name> <nr of cells to plot> filename"
  sys.exit(1)
else:
  figname=sys.argv[1]
  tracknr=int(sys.argv[2])

filecounter=0
for filecounter,filename in enumerate(sys.argv[3:]):

  print "reading file ",filename
  timepoints=[]
  speed=[]
  avspeed=[]
  stdevspeed=[]
  xpos=[]
  ypos=[]
  xpos.append([])
  ypos.append([])
  autocorr=[]
  
  count=0
  
  ##read data from file: store raw cell positions
  with open(filename,"r") as fin:

    #read first line first to get the first time point (could probably more clever, but hey)
    for line in fin:
      line=line.split()
      timepoints.append(int(line[0]))
      xpos[-1].append(float(line[3]))
      ypos[-1].append(float(line[4]))
      count+=1
      break

    #read rest of file
    for line in fin:
      line=line.split()
      if int(line[0])!=timepoints[-1]:
        xpos.append([-1.]*len(xpos[0]))
        ypos.append([-1.]*len(xpos[0]))
        count=0
        timepoints.append(int(line[0]))
      if (timepoints[-1]==timepoints[0]):
        xpos[-1].append(float(line[3]))
        ypos[-1].append(float(line[4]))
      else:
        xpos[-1][int(line[1])-1]=float(line[3])
        ypos[-1][int(line[1])-1]=float(line[4])
      count+=1

  timeint=timepoints[-1]-timepoints[-2]
  lminy=[]
  lavy=[]
  maxint=len(xpos)
  for i in range(1,maxint):
    miny = min(ypos[i])
    lminy.append(miny)
    avy=np.mean(ypos[i])
    lavy.append(avy)
    # if i>1: 
    #     # print avxpos
    #     squarex = (lcmx[-1] - lcmx[-2])**2
    #     squarey = (lcmy[-1] - lcmy[-2])**2
    #     dist = math.sqrt(  squarex  +  squarey  )
    #     lcmspeed.append( ( 1./float(timeint) ) * dist  )

  ##start plotting##
  #average speed through simulation
  #ax0.plot(timepoints,MSD)
  #print stdevspeed
  print "last time point:",timepoints[-1]   
  #label = filename.split('.')[-2].split('_')[-2:-1]
  ax0.plot(timepoints[1:],lminy,c=colours[filecounter],lw=0.5,label=filename)
  #ax0.set_yscale('log')
  ax0.legend(fontsize='xx-small',title='first cell')
  ax1.plot(timepoints[1:],lavy,c=colours[filecounter],lw=0.5)
  #ax1.legend(title='center of mass')

#print "axis: ", ax0.get_xaxis()
#ax0.set_xlim([0, 150000])
#ax1.set_xlim([0, 150000])
plt.xlim([0,150000])
fig.savefig(figname, bbox_inches='tight')
# plt.show()
