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

if len(sys.argv) <2:
    print ("This is the program 'plot_division_time.py'")
    print ("Usage: ./plot_division_time.py <figure name> <season duration> filename")
    sys.exit(1)
else:
    figname=sys.argv[1]
    season=int(sys.argv[2])
    filename=sys.argv[3]

##read data from file: store division data
with open(filename,"r") as fin:

    #read first line first to get the first time point (could probably more clever, but hey)
    divtime=[]
    divtime.extend([] for i in range(3))
    medians=[]
    medians.extend([] for i in range(3))
    seasoncount=1
    far=0
    add=0
    cells=[0 for x in range(5000)]
    for line in fin:
      line=line.split(' ')
      prevtime=int(line[0])
      cells[int(line[1])]=int(line[11])
      break

    print ("reading rest of file")
    #read rest of file
    for line in fin:

        line=line.split(' ')
        time=int(line[0])
        thiscell=int(line[1])
        divnum=int(line[11]) 

        #new season: plot previous and reset vars
        if time!=prevtime:
            prevtime=time
            if not(prevtime%season):
                if divnum!=0:
                    print ("oh oh...{} {}".format(divnum, prevtime))
                    sys.exit(1)
                for ind,div in enumerate(divtime):        
                    #print ind, np.median(div)
                    if (len(div)):      
                        medians[ind].append(np.median(div))
                    else: 
                        print ("no median for div {}".format(ind))
                        medians[ind].append(0)
                if(far):
                    ax0.scatter(prevtime/season, 0, c="black", s=5)
            
                #reset vars
                divtime=[]
                divtime.extend([] for i in range(3))
                cells=[0 for x in range(5000)]
                seasoncount+=1
                far=0
                
        #check if the cells have to travel far
        if(not far and int(line[12])<300):
            far=1

        if divnum>cells[thiscell]:
            #print divnum
            divtime[divnum-1].append(time%season)
            cells[thiscell]=divnum
           

for ind, div in enumerate(medians):
    #print ind
    ax0.plot(range(1,len(div)+1), div, color=colours[ind])
ax0.set_xlabel('season')
ax0.set_ylabel('division time')
ax0.set_title('timing of cell divisions in a season')

fig.savefig(figname, bbox_inches='tight')
#plt.show()
