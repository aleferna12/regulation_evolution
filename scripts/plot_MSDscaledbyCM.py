#!/usr/bin/python2.7

'''
Plot the MSD over time period, minus CM for one or more simulations.
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
  print "This is the program 'plot_MSD_blabla.py'"
  print "Usage: ./plot_MSD_blabla_.py <figure name> <nr of tracks to plot> <filename(s)> "
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
  xpos=[] # [ [list of xpos at first timestep], [list of xpos at second timestep], ... ]
  ypos=[]
  xn=[]
  xpos.append([])
  ypos.append([])
  count=0
  
  xminCM=[]
  yminCM=[]
  MSDminCM=[]
  SDEVminCM=[]
  
  output = subprocess.Popen(['head', '-1', filename], stdout=subprocess.PIPE).communicate()[0]
  timepoints.append( int(output.split(' ')[0]) )
  
  ##read data from file: store raw cell positions
  with open(filename,"r") as fin:
    #read file
    for line in fin:
      
      line=line.split()
      if int(line[0]) != timepoints[-1]:
        #get xpos -> calculate CM, subtract from xpos
        xCM = np.mean(xpos[-1])
        yCM = np.mean(ypos[-1])
        xminCM.append([ x - xCM for x in xpos[-1] ])
        yminCM.append([ y - yCM for y in ypos[-1] ])
        
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
  
  print "length xpos", len(xpos)
  print "length xminCM", len(xminCM)
  
  ##calculate MSD and standard error
  maxint=len(xpos) -1 #because we don't do calculations on the last time point of xminCM and yminCM
  nrcells=len(xpos[0])
  print maxint
  factor_for_data = 100 # set this to 1 to get all the data
  for a in range(maxint/factor_for_data):
    i = a*factor_for_data #do it every 10 time steps
    print i
    if i> 20000/(timepoints[-1] - timepoints[-2]) : break
    
    count=0
    count2=0
    MSD.append(0.0)
    MSDminCM.append(0.0)
    sd=0.0
    sd2=0.
    xn=[]
    xn2=[]
    deltaxminCM=0
    deltayminCM=0
    count_count=0
    for c in range(nrcells):  #problem when cells die...
        deltaxminCM += xminCM[count+1][c]-xminCM[count][c]
        deltayminCM += yminCM[count+1][c]-yminCM[count][c]
        count_count +=1
    
    while (count+i< maxint):
      
      for c in range(nrcells):  #problem when cells die...
        if xpos[count+i][c]>-1:
          xn.append((xpos[count+i][c]-xpos[count][c])**2 + (ypos[count+i][c]-ypos[count][c])**2 )
          xn2.append( (xminCM[count+i][c]-xminCM[count][c])**2 + (yminCM[count+i][c]-yminCM[count][c])**2 )
          
          MSD[-1]+=xn[-1]
          MSDminCM[-1]+=xn2[-1]
          count2+=1
        else:
          print "cell ",c,"is no more"
      count+=1
    
    MSD[-1]=MSD[-1]/float(count2)
    MSDminCM[-1]=MSDminCM[-1]/float(count2)
    for el in xn:
      sd+=(el-MSD[-1])*(el-MSD[-1])
    for el in xn2:
      sd2+=(el-MSDminCM[-1])**2.
    SDEV.append(math.sqrt(sd/float(count2)))
    SDEVminCM.append(math.sqrt(sd2/float(count2)))

  #center cell track start at 0,0  
  # print "cell tracks"
  # for i in range(1,len(xpos)):
  #   for j in range(len(xpos[0])):
  #     if xpos[i][j]>-1: 
  #       xpos[i][j]-=xpos[0][j]
  #       ypos[i][j]-=ypos[0][j] 
  #     else:
  #       xpos[i][j]=xpos[i-1][j]
  #       ypos[i][j]=ypos[i-1][j] 
  # 
  # for j in range(len(xpos[0])): 
  #   xpos[0][j]-=xpos[0][j]
  #   ypos[0][j]-=ypos[0][j] 
  
  ##start plotting##
  
  
  #ax0.plot(timepoints,MSD)
  print "len timepoitns and MSD = ", len(timepoints[:-1:factor_for_data]), len(MSD)
  
  print "average deltaxminCM", deltaxminCM/float(count_count)
  print "average deltayminCM", deltayminCM/float(count_count)
  
  deltayminCM += yminCM[count+i][c]-yminCM[count][c]
  count_count +=1
  
  # print timepoints[:-1:factor_for_data]
  # print MSD
  ax0.errorbar(timepoints[:-1:factor_for_data][:len(MSD)],MSD, yerr=SDEV, fmt='', c=colours[filecounter],errorevery=100)
  ax0.errorbar(timepoints[:-1:factor_for_data][:len(MSD)],MSDminCM, yerr=SDEVminCM, fmt='', c=colours[filecounter+1],errorevery=100)
  #ax0.set_yscale('log')


  #plot cell tracks
  #invert position matrices for plotting (each row is 1 cell instead of 1 timepoint)
  xx=zip(*xpos)
  yy=zip(*ypos)
  toplot=range(0,len(xx))
  random.shuffle(toplot)
  #for el in toplot[:tracknr]:
  #  print el

  count=0
  count2=0
  for xrow,yrow in zip(xx,yy):
    if count in toplot[:tracknr]:
      ax1.plot(xrow,yrow, lw=0.5 ,c=colours[filecounter])
      count2+=1
    count+=1

  filecounter+=1

ax0.set_xlabel('time (MCS)')
ax0.set_ylabel('MSD (pix^2)')
ax0.set_title('Mean squared displacement')

ax1.set_xlabel('x position')
ax1.set_ylabel('y position')
ax1.set_title('Re-centered cell tracks')

fig.savefig(figname, bbox_inches='tight')
plt.show()
