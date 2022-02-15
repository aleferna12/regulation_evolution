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

if len(sys.argv) <2:
    print ("This is the program 'plot_division_time_gatherseasons.py'")
    print ("Usage: ./plot_division_time.py <figure name> <which division (1, 2, 3)> filename")
    sys.exit(1)
else:
    figname=sys.argv[1]
    div=int(sys.argv[2])
    filename=sys.argv[3]

##read data from file: store division data
with open(filename,"r") as fin:

    #read first line first to get the first time point (could probably more clever, but hey)
    labels=[]
    xpos=[]
   
    for line in fin:
      line=line.split(' ')
      labels.append(line[0])
      xpos.append(float(line[div]))
      


ax0.set_xlabel('average moment of first division')
ax0.set_ylabel('division time')
ax0.set_title(os.path.splitext(filename)[0])
ax0.scatter(xpos, [1]*len(xpos), marker="o", s=50)
for ind,el in enumerate(labels):
    ypos=1+random.uniform(-0.001, 0.001)
    print ("{}".format(ypos))
    ax0.annotate(el, (xpos[ind]-500,ypos))

fig.savefig(figname, bbox_inches='tight')
#plt.show()
