#!/usr/bin/python2.7

'''
Plot the MSD over time for one or more simulations.
Plot some selected cell tracks, centered at 0,0
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
fig, (ax0, ax1) = plt.subplots(ncols=2)
#fig, ax0 = plt.subplots()

if len(sys.argv) <3:
  print "This is the program 'plot_MSD2.py'"
  print "Usage: ./plot_MSD.py <figure name> <nr of tracks to plot> <filename(s)> "
  sys.exit(1)
else:
  figname=sys.argv[1]
  tracknr=int(sys.argv[2])

filecounter=0
for filename in sys.argv[3:]:

  print "reading file ",filename
  timepoints=[]
  MSD=[]
  SDEV=[]
  xpos=[]
  ypos=[]
  xn=[]
  xpos.append([])
  ypos.append([])
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
  
  #I should eliminate the initial data, because it is initial condition
  toremove=len(timepoints)/4
  xpos = xpos[toremove:]
  ypos = ypos[toremove:]
  timepoints=timepoints[: -toremove ] #this has to be like this because timepoints is also used as x axis
  
  print "Done reading, now calculating MSD"
  ##calculate MSD and standard error
  maxint=len(xpos)
  nrcells=len(xpos[0])
  print maxint
  for i in range(maxint):
    count=0
    count2=0
    MSD.append(0.0)
    sd=0.0
    xn=[]
    while (count+i<maxint):
      
      for c in range(nrcells):  #problem when cells die...
        if xpos[count+i][c]>-1:
          xn.append((xpos[count+i][c]-xpos[count][c])*(xpos[count+i][c]-xpos[count][c])+(ypos[count+i][c]-ypos[count][c])*(ypos[count+i][c]-ypos[count][c]))
          MSD[-1]+=xn[-1]
          count2+=1
        else:
          print "cell ",c,"is no more"
      count+=1
    MSD[-1]=MSD[-1]/float(count2)
    for el in xn:
      sd+=(el-MSD[-1])*(el-MSD[-1])
    SDEV.append(math.sqrt(sd/float(count2)))

  #center cell track start at 0,0  
  print "Done that, now cell tracks"
  for i in range(1,len(xpos)):
    for j in range(len(xpos[0])):
      if xpos[i][j]>-1: 
        xpos[i][j]-=xpos[0][j]
        ypos[i][j]-=ypos[0][j] 
      else:
        xpos[i][j]=xpos[i-1][j]
        ypos[i][j]=ypos[i-1][j] 

  for j in range(len(xpos[0])): 
    xpos[0][j]-=xpos[0][j]
    ypos[0][j]-=ypos[0][j] 

  ##start plotting##

  print "timepoints: ", timepoints
  print "MSD:", MSD
  
  # ax0.set_xscale('log')
  # ax0.set_yscale('log')
  #to linear interpolation of log log transformed data, 
  # to get exponent
  log_t=timepoints[:]
  log_MSD=MSD[:]
  if log_t[0]==0:
      log_t=log_t[1:]
      log_MSD=log_MSD[1:]
  
  log_t = [np.log(x) for x in log_t]
  log_MSD = [np.log(x) for x in log_MSD]
  m_and_b = np.polyfit(log_t, log_MSD, 1)
  poly1d_fn = np.poly1d(m_and_b) # CAREFUL HERE: poly1d_fn is now a function which takes in x and returns an estimate for y
  print "filename: ", filename, "m,b:", m_and_b
  ax0.plot(timepoints,MSD)
  ax1.plot(log_t,log_MSD, 'yo', markersize=0.5) 
  ax1.plot(log_t, poly1d_fn(log_t), '--k')
  
  ax0.errorbar(timepoints,MSD, yerr=SDEV, fmt='', c=colours[filecounter],errorevery=100)
  
  #ax0.set_yscale('log')


  # #plot cell tracks
  # #invert position matrices for plotting (each row is 1 cell instead of 1 timepoint)
  # xx=zip(*xpos)
  # yy=zip(*ypos)
  # toplot=range(0,len(xx))
  # random.shuffle(toplot)
  # #for el in toplot[:tracknr]:
  # #  print el
  # 
  # count=0
  # count2=0
  # for xrow,yrow in zip(xx,yy):
  #   if count in toplot[:tracknr]:
  #     ax1.plot(xrow,yrow, lw=0.5 ,c=colours[filecounter])
  #     count2+=1
  #   count+=1

  filecounter+=1

ax0.set_xlabel('time (MCS)')
ax0.set_ylabel('MSD (pix^2)')
ax0.set_title('Mean squared displacement')

ax1.set_xlabel('log(time)')
ax1.set_ylabel('log(MSD)')
ax1.set_title('log-log-plot(MSD/time)')

fig.savefig(figname, bbox_inches='tight')
#plt.show()
