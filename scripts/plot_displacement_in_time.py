#!/usr/bin/python2.7

'''
Plot position in time, takes all cells in a time step and plots their average position
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

colours=[(0.7,0.13,0.),"royalblue", "darkgoldenrod", (0.,0.5,0.2), "salmon", "lightskyblue","orchid"]
filename=""
fig, (ax0, ax1) = plt.subplots(ncols=2,sharey=True)
#fig, ax0 = plt.subplots()
# av_xpos=[[],[],[]]
# av_ypos=[[],[],[]]
# lcount=[[],[],[]]

if len(sys.argv) <3:
  print "This is the program 'plot_diplacement_in_time'"
  print "Usage: ./plot_speed.py <figure name> filename"
  sys.exit(1)
else:
  figname=sys.argv[1]
  #tracknr=int(sys.argv[2])

peakx = 500
peaky = 0
print "Calculating distance from peak at (row,col) = ", peakx,peaky

filecounter=0
for filename in sys.argv[2:]:

    print "reading file ",filename
    # timepoints=[]
    # speed=[]
    # avspeed=[]
    # stdevspeed=[]
    # xpos=[]
    # ypos=[]
    # xpos.append([])
    # ypos.append([])
    # autocorr=[]

    # count=0

    av_xpos=[]
    av_ypos=[]
    min_xpos=[]
    min_ypos=[]
    ltime=[]
    lxpos=[]
    lypos=[]
    c1_pos=[]
    ldist=[]
    av_dist=[]
    min_dist=[]
    
    # read time in the first line
    output = subprocess.Popen(['head', '-1', filename], stdout=subprocess.PIPE).communicate()[0]
    time = int(output.split(' ')[0])
    ltime = [time]        
    with open(filename,"r") as fin:
        #read first line first to get the first time point (could probably more clever, but hey)
        for line in fin:
            line=line.split()
      
            time = int(line[0])
            xpos = float(line[3])
            ypos = float(line[4])
            dist = np.hypot( xpos-peakx,ypos-peaky )
            # print xpos,ypos
            #if time of this line different from time of previous line
            # we do aggreagate statistics
            if time != ltime[-1]:
                
                av_xpos.append(  np.mean(  lxpos )  )
                min_xpos.append( min(lxpos) )
                lxpos=[]
                
                av_ypos.append(  np.mean(  lypos )  )
                min_ypos.append( min(lypos) )
                lypos=[] 
                
                av_dist.append( np.mean(ldist) )
                min_dist.append( min(ldist))
                ldist=[]
                
                #print ltime[-1], min_xpos[-1],av_xpos[-1]
                ltime.append(time)
              
            lxpos.append(xpos)
            lypos.append(ypos)
            ldist.append(dist)
          
    printthis=0
    
    #colouring and data handling
    # if 'CE_g-4' in filename: 
    #     index=0
    #     printthis=1
    # elif 'CE_g0' in filename: 
    #     index=1
    #     printthis=1
    # elif 'CE_g4' in filename: 
    #     index=2
    #     printthis=1
    # elif 'CE_g6' in filename: 
    #     index=3
    #     printthis=1
    # else: 
    #     index=4
    #     printthis = 0 # zeeeeeroooooo
    # 
    if 'g6_1c_f1000' in filename: 
        index=0
        printthis=1
    elif 'g6_50c_f1000' in filename:
        index=1
        printthis=1
    else:
        index=-1
        printthis=1
    # if 1:
    if printthis:
        #disregard last time point because algotirthm above exsits by then
        
        # ax0.plot( ltime[:-1],av_ypos, c= colours[index], lw=0.5  )
        # ax1.plot( ltime[:-1],min_ypos, c= colours[index], lw=0.25  )
        # 
        ax0.plot( ltime[:-1],av_dist, c= colours[index], lw=0.5  )
        ax1.plot( ltime[:-1],min_dist, c= colours[index], lw=0.25  )
    
    
    if not printthis:
        ax1.plot( ltime[:-1],min_ypos, c= colours[index], lw=0.25  )
        if not c1_pos:
            c1_pos = av_ypos
            howmany1c=1
            ltime_1c = ltime[:]
        else:
            c1_pos = [ x+y for x,y in zip(c1_pos, av_ypos) ]   #zp works only over the smallest list, which is exactly what we want
            howmany1c+=1
            if len(ltime)<ltime_1c:
                ltime_1c = ltime[:]
            
try:
  ax0.plot( ltime_1c[:-1],[x/float(howmany1c) for x in c1_pos], c= colours[4], lw=0.5  )
except:
    pass

  # ax0.plot( ltime,av_xpos, c=colour, lw=0.2  )
  # ax1.plot( ltime,min_xpos, c=colour, lw=0.2  )

# print [len(x) for x in av_xpos]    
# for i in range(len(av_xpos)):
#     print i
#     ax0.plot( timepoints[:len(av_xpos[i])],[ x/float(lcount[i]) for x in av_xpos[i] ], c=colours[i], lw=1.  )
#     ax1.plot( timepoints[:len(av_ypos[i])],[ y/float(lcount[i]) for y in av_ypos[i] ], c=colours[i], lw=1.  )


# ax0.set_xlabel('time (MCS)')
# ax0.set_ylabel('average cell speed (pix/MCS)')
# ax0.set_title('Average speed through simulation')
# if not b_centerofmass:
#     ax0.set_xlim(0., 20000.)
# 
# ax1.set_xlabel('lag (MCS)')
# ax1.set_ylabel('autocorrelation')
# #ax1.set_ylabel('instantaneous speed (pix/MCS)')
# ax1.set_title('NOT Speed autocorrelation')
# ax1.set_xlim(0., 20000.)
# 

ax0.set_xlim([0,150000])
ax0.set_ylim([0,None])
ax1.set_xlim([0,150000])
ax1.set_ylim([0,None])
fig.savefig(figname, bbox_inches='tight')
# plt.show()
