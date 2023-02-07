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

colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
fig, (ax0, ax1) = plt.subplots(nrows=2)
#fig, ax0 = plt.subplots()

if len(sys.argv) <3:
  print "This is the program 'plot_speed.py'"
  print "Usage: ./plot_speed.py <figure name> <nr of cells to plot> filename"
  sys.exit(1)
else:
  figname=sys.argv[1]
  tracknr=int(sys.argv[2])

filecounter=0
for filename in sys.argv[3:]:

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

  ##calculate MSD and standard error
  timeint=timepoints[-1]-timepoints[-2]
  maxint=len(xpos)
  nrcells=len(xpos[0])
  print maxint, nrcells
  count2=0
  totav=0.
  cellstdev=[0.0]*nrcells
  
  b_centerofmass = True
  lcmx = []
  lcmy = []
  lcmspeed =[]
  
  if not b_centerofmass:
      for i in range(1,maxint):
        speed.append([])
        avspeed.append(0.0)
        sd=0.0
        count=0
        for c in range(nrcells):  #problem when cells die...
          if xpos[i][c]>-1:
            displace=math.sqrt( (xpos[i][c]-xpos[i-1][c])*(xpos[i][c]-xpos[i-1][c])+(ypos[i][c]-ypos[i-1][c])*(ypos[i][c]-ypos[i-1][c]) ) / float(timeint)
            speed[-1].append(displace)
            avspeed[-1]+=displace
            totav+=displace
            count+=1
            count2+=1
          else:
            print "cell ",c,"is no more"
            speed[-1].append(-1)
        avspeed[-1]=avspeed[-1]/float(count)
        cc=0
        for el in speed[-1]:
          if el>-1:
            sd+=(el-avspeed[-1])*(el-avspeed[-1])
            cellstdev[cc]+=(el-avspeed[-1])*(el-avspeed[-1])
          cc+=1
        stdevspeed.append(math.sqrt(sd/float(count)))
  else:
        for i in range(1,maxint):
            avxpos = np.mean(xpos[i])
            avypos = np.mean(ypos[i])
            lcmx.append(avxpos)
            lcmy.append(avypos)
            if i>1: 
                # print avxpos
                squarex = (lcmx[-1] - lcmx[-2])**2
                squarey = (lcmy[-1] - lcmy[-2])**2
                dist = math.sqrt(  squarex  +  squarey  )
                lcmspeed.append( ( 1./float(timeint) ) * dist  )

  #overall average
  if not b_centerofmass:
      totav/=float(count2)
      print "average speed = ",totav 

      for i in range(maxint-1): #lag period
        ccc=0
        autocorr.append([0.0]*nrcells)
        while (ccc+i<maxint-1):
          count2=0
          for c in range(nrcells):  
            if speed[ccc+i][c]>-1:
              autocorr[-1][c]+=(speed[ccc][c]-totav)*(speed[ccc+i][c]-totav)
          ccc+=1  
        for c in range(nrcells):  
          autocorr[-1][c]/=cellstdev[c]
      

  ##start plotting##
  #average speed through simulation
  #ax0.plot(timepoints,MSD)
  #print stdevspeed
  if not b_centerofmass:
    ax0.errorbar(timepoints[1:],avspeed, yerr=stdevspeed, fmt='-o', c=colours[filecounter],errorevery=20)
  else:
    print "last time point:",timepoints[-1]   
    ax0.plot(timepoints[2:],lcmspeed,c=colours[filecounter],lw=0.5)
  #ax0.set_yscale('log')
  

  #plot individual cell's speeds
  #invert speed matrix for plotting (each row is 1 cell instead of 1 timepoint)
  cellspeed=zip(*speed)
  toplot=range(0,len(cellspeed))
  random.shuffle(toplot)

  ##autocorrelation plotting
  #cellauto=zip(*autocorr)
  #toplot=range(0,len(cellauto))
  #random.shuffle(toplot)

  count=0
  for cel in cellspeed:
    if count in toplot[:tracknr]:
      ax1.plot(timepoints[1:],cel, c=colours[filecounter], alpha=0.5)
    count+=1

  ##speed histograms
  #flat_speed = []
  #flat_speed = [y for x in speed for y in x]
  #ax1.hist( flat_speed,50,normed=1 )

  
  filecounter+=1


ax0.set_xlabel('time (MCS)')
ax0.set_ylabel('average cell speed (pix/MCS)')
ax0.set_title('Average speed through simulation')
if not b_centerofmass:
    ax0.set_xlim(0., 20000.)

ax1.set_xlabel('lag (MCS)')
ax1.set_ylabel('autocorrelation')
#ax1.set_ylabel('instantaneous speed (pix/MCS)')
ax1.set_title('NOT Speed autocorrelation')
ax1.set_xlim(0., 20000.)



fig.savefig(figname, bbox_inches='tight')
# plt.show()
