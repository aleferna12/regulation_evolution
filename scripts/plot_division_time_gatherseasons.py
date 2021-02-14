#!/usr/bin/python3.8

'''
plot the times at which cells divide in the simulations
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
#fig, (ax0, ax1) = plt.subplots(nrows=2)
fig, ax0 = plt.subplots()

if len(sys.argv) <5:
    print ("This is the program 'plot_division_time_gatherseasons.py'")
    print ("Usage: ./plot_division_time.py <figure name> <season duration> <startseason> <endseason> filename")
    sys.exit(1)
else:
    figname=sys.argv[1]
    season=int(sys.argv[2])
    start=int(sys.argv[3])
    end=int(sys.argv[4])
    filename=sys.argv[5]

##read data from file: store division data
with open(filename,"r") as fin:

    #read first line first to get the first time point (could probably more clever, but hey)
    divtime=[]
    divtime.extend([] for i in range(3))
    seasoncount=1
    far=0
    add=0
    cells=[0 for x in range(5000)]
    for line in fin:
      line=line.split(' ')
      prevtime=int(line[0])
      if (prevtime==start): 
          cells[int(line[1])]=int(line[11])
          break

    # print ("reading the rest")
    #read rest of file
    for line in fin:

        line=line.split(' ')
        time=int(line[0])
        thiscell=int(line[1])
        divnum=int(line[11]) 

        #new season: plot previous and reset vars
        if time==end:
            break
        if time!=prevtime:
            prevtime=time
            if not(prevtime%season):
                if divnum!=0:
                    print ("oh oh...{}, {}".format(divnum, prevtime))
                    sys.exit(1)
        
                #if(far):
                #    print ("it's far")
            
                #reset vars
                #divtime=[]
                #divtime.extend([] for i in range(3))
                cells=[0 for x in range(5000)]
                far=0
                seasoncount+=1
                
        #check if the cells have to travel far
        if(not far and int(line[12])<300):
            far=1

        if divnum>cells[thiscell]:
            #print divnum
            divtime[divnum-1].append(time%season)
            cells[thiscell]=divnum
           
#print ("nr seasons is {}".format(seasoncount))
print ("{} {}".format(os.path.splitext(filename)[0],sum(divtime[0])/float(len(divtime[0])) ) )
div1=ax0.violinplot(divtime,positions=[1,2,3],widths=0.5, showmedians=True, showextrema=False)
#div1['bodies'].set_facecolor(colours[:3])
div1['cmedians'].set_color([ colors.to_rgba(x) for x in colours[:3]  ]  )
for ind,pc in enumerate(div1['bodies']):
    pc.set_facecolor(colours[ind])

ax0.set_ylim(0,500000)
ax0.set_xlabel('division number')
ax0.set_ylabel('division time')
ax0.set_title(os.path.splitext(filename)[0])

fig.savefig(figname, bbox_inches='tight')
#plt.show()
