#!/usr/bin/python2.7

'''
Plot the average angle between persistence vectors of cells given a particular distance between these cells
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
###     function      ###
###                   ###
#########################

def blob_plot(ax,colour,data,list_pos, extra=False):
    lextra=[]
    for d,p in zip(data,list_pos):
      
      #check if there are data in the lists
      if len(d):
        #this makes the histogram that constitutes the blob plot
        m=min(d)
        M=max(d)
        nbins=30
        x = np.linspace(m,M,nbins) # support for histogram
        his,bins = np.histogram(d, bins=nbins)
      
        #scale the histogram to a proper size (you may have to fiddle with this), then place it at position 'pos'
        scale_blob=3.
        shift_his_plus_pos =  [ p + h*scale_blob/float(len(d))  for h in his]
        shift_his_minus_pos = [ p - h*scale_blob/float(len(d))  for h in his]
      
        facecolor,alpha=(colour,0.5) # assign color and transparency
        #this is the matplotlib function that does the trick
        ax.fill_betweenx(x,shift_his_minus_pos, shift_his_plus_pos, linewidth=0.0, facecolor= facecolor, alpha=alpha)
        #calculates median or mean, if you want
        if extra=='median': lextra.append( np.median(d) )
        elif extra=='mean': lextra.append( np.mean(d) )
        #and plots it
        if extra != False:
          color = 'orangered'
          ax.plot(list_pos,lextra, color=color_mean_or_median,linestyle='-',marker='D',markersize=5, lw=1.5)
      


#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
fig, (ax0,ax1) = plt.subplots(nrows=2,sharex=True)
fig.subplots_adjust(hspace=0)

#fig, ax0 = plt.subplots()

if len(sys.argv) <3:
  print "This is the program 'plot_angle_correlations.py'"
  print "Usage: ./plot_angle_correlations.py <figure name> <cell targetarea> <binsize> <blobplots> <filename(s)> "
  sys.exit(1)
else:
  figname=sys.argv[1]
  celldiam=math.sqrt(float(sys.argv[2])/math.pi)*2.
  print "cell diameter is ",celldiam
  binsize=float(sys.argv[3])
  binmult=1./binsize
  blobs=int(sys.argv[4])

maxdist=0
filecounter=0


for filename in sys.argv[5:]:

  print filename
  avangle=[0.]*200*int(binmult)   #1 bin is binsize celldiam 
  binamount=[0]*200*int(binmult)  #to keep track of how many in each bin
  angles=[[] for x in range(200*int(binmult))]
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
      time0=int(line[0])
      timepoint=time0
      xpos[-1].append(float(line[3]))
      ypos[-1].append(float(line[4]))
      count+=1
      break

    #read rest of file
    for line in fin:
      line=line.split()
      if int(line[0])!=timepoint:
        xpos.append([-1.]*len(xpos[0]))
        ypos.append([-1.]*len(xpos[0]))
        count=0
        timepoint=int(line[0])
      if timepoint==time0:
        xpos[-1].append(float(line[3]))
        ypos[-1].append(float(line[4]))
      else:
        xpos[-1][int(line[1])-1]=float(line[3])
        ypos[-1][int(line[1])-1]=float(line[4])
      count+=1

  ##calculate MSD and standard error
  maxint=len(xpos)
  nrcells=len(xpos[0])
  print "maxint=",maxint
  for i in range(500,maxint): #skip first time step, which is weird
    count=0
    count2=0
    xn=[]      
    for c1 in range(nrcells):  #problem when cells die...
      for c2 in range(nrcells):
        if c1 != c2 and xpos[i][c1]>-1 and xpos[i][c2]>-1:
          dist=math.sqrt((xpos[i][c1]-xpos[i][c2])**2+(ypos[i][c1]-ypos[i][c2])**2)
          dist2=math.sqrt((xpos[i-1][c1]-xpos[i-1][c2])**2+(ypos[i-1][c1]-ypos[i-1][c2])**2)
          dist=(dist+dist2)/2.
          #print "dist is ", dist, "bin nr is ",int(dist/celldiam*100.)
          xmov1=xpos[i][c1]-xpos[i-1][c1]
          ymov1=ypos[i][c1]-ypos[i-1][c1]
          xmov2=xpos[i][c2]-xpos[i-1][c2]
          ymov2=ypos[i][c2]-ypos[i-1][c2]
          if ( math.sqrt(xmov1**2+ymov1**2) * math.sqrt(xmov2**2+ymov2**2) )>0.00000:
            dotprod=(xmov1*xmov2+ymov1*ymov2) / ( math.sqrt(xmov1**2+ymov1**2) * math.sqrt(xmov2**2+ymov2**2) )
          else:
            dotprod=-2.
          if dotprod<1. and dotprod>-1.:
            avangle[int(dist/celldiam*binmult+binmult*0.5)]+=np.arccos(dotprod)*180/math.pi
            angles[int(dist/celldiam*binmult+binmult*0.5)].append(np.arccos(dotprod)*180/math.pi)
            binamount[int(dist/celldiam*binmult+binmult*0.5)]+=1
          else:
            print "problem: dot product is ","{0:.5f}".format(dotprod)
  xpoints=[]
  start=-1
  cont=0
  print avangle[:15]
  totbins=0
  for ind, el in enumerate(binamount):
    if el>0:
      if start==-1:
        start=ind
      avangle[ind]/=float(el)
      xpoints.append(float(ind)*binsize)
      totbins+=el
      cont=0
    elif start!=-1 and el==0:
      cont+=1
      avangle[ind]=-1
      xpoints.append(float(ind)*binsize)
      if cont>5:
        final=ind
        if final>maxdist:
          maxdist=final
        break
  
  #print start, final
  #print xpoints
  
  binnorm=[float(x)/float(totbins) for x in binamount]

  if blobs!=0:
    blob_plot(ax0, colours[filecounter], angles[start:], list_pos=xpoints)	#makes blob plots, if extra is not spcified, makes no median/mean
  
  
  #colouring and data handling
  if 'CE_g-4' in filename: 
      index=0
      printthis=1
  elif 'CE_g0' in filename: 
      index=1
      printthis=1
  elif 'CE_g4' in filename: 
      index=2
      printthis=1
  elif 'CE_g6' in filename: 
      index=3
      printthis=1
  else: 
      index=4
      #printthis = 0 # zeeeeeroooooo
  
  ax0.plot(xpoints[:-5],avangle[start:final-4], linewidth=0.75, color=colours[index])
  ax1.plot(xpoints[:-5],binamount[start:final-4],marker='o', linewidth=0., ms=0.6, color=colours[index],alpha=0.75)

  filecounter+=1

x = np.linspace(0,int(maxdist*binsize),100)
y = 0*x+90
ax0.plot(x,y,linestyle='dashed',color="black", linewidth=0.5)
ax0.set_xlim(0,8)
ax0.set_ylim(45,120)
ax1.set_xlabel('distance (cell diameters)')
ax0.set_ylabel('angle (degrees)')
ax0.set_title('angle between displacement directions of cell pairs')
ax1.set_ylabel('fraction of observations')

fig.savefig(figname, bbox_inches='tight')
#plt.show()
