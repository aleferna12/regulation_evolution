#!/usr/bin/python3.8

'''
plot the times at which cells divide in the simulations
'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np
import seaborn as sns
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

if len(sys.argv) <2:
    print ("This is the program 'plot_divtimechange_allsims.py'")
    print ("Usage: ./plot_divtimechange_allsims.py <figure name> <simnum> <filenames....>")
    sys.exit(1)
else:
    figname=sys.argv[1]
    simnum=int(sys.argv[2])

    
#initialise figure
fig, ax = plt.subplots()
colorset=sns.color_palette("twilight",n_colors=simnum)
pos1 = 0
pos2 = 2
pos3 = 4

startdivs={} #what was the average division time at the start of the simulation
startpos={} #where to place the initial dot
colcount=0
with open(sys.argv[3]) as fi:
    for line in fi:
        line=line.split(' ')
        div= float(line[1])  #float(line[2]) line 1 is av, line2 is mdeian
        check=0
        for el in startdivs.values():
            if math.sqrt((div-el)**2)<1000:
                check+=1
        startpos[int(line[0])]=pos2-0.05+check*0.05
        startdivs[int(line[0])]=div
        
        

##read data from file: store division data
count=0
colcount=-1
previouswhich=-1
for filename in sys.argv[4:]:
    if (not(count%3)):
        whichstart=int(filename)
        if previouswhich is not whichstart:
            colcount+=1
            colcount=colcount%simnum
            previouswhich=whichstart
    elif (count%3)==1:
        pos=int(filename) #whether it is evolved with adhesion or without
        #print (pos)
    else:
        with open(filename,"r") as fin:

            #data should be straightforward
            for line in fin:
                arr=line.split()
                enddiv=float(arr[1])
                secpos=pos*4
                #print (whichstart, colcount, pos)
                ax.plot( [ startpos[whichstart], secpos+0.2*(random.random()-0.5) ] , [startdivs[whichstart],enddiv], \
                lw=1, alpha=0.7, marker='.', markersize = 15, color=colorset[colcount])

    count+=1

#ax.set_ylim(0,250000)

x = np.array([0,2,4])
my_xticks = ['without\nadhesion','ancestor','with\nadhesion']
plt.xticks(x, my_xticks)

#ax0.set_xlabel('evolution')
ax.set_ylabel('first division time')
fig.savefig(figname, bbox_inches='tight')
plt.show()
