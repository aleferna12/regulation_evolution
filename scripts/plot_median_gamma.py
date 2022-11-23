#!/usr/bin/python3.8

'''
plot the average gamma value over time
'''

import matplotlib as mpl
mpl.use('QT5Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import sys,math,os,subprocess,random
import seaborn as sns
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
#from PIL import Image


#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
fig, ax0 = plt.subplots()


## Check commandline arguments ##
if len(sys.argv) <2:
    print ("This is the program 'plot_median_gamma.py'")
    print ("Usage: ./plot_median_gamma.py <figure name> <season duration>")
    sys.exit(1)
else:
    figname=sys.argv[1]
    season=int(sys.argv[2])

nrfiles=len(sys.argv)-3

## make a color palette to color the trends of different simulations (seaborn has more pretty palettes)
gamcolorset=sns.color_palette("Blues",n_colors=nrfiles)


#########################
## Begin reading files ##
#########################

filecount=0
for filename in sys.argv[3:]:
    print(filename)
    filecount+=1
    
    ## read data from file: store adhesion data ##
    with open(filename,"r") as fin:

        #initialise storage lists
        adval=[]
        adval.extend([] for i in range(3)) #in this case: median, low, high
        storetime=[]
        seasoncount=1
        add=0
        lJmed=[]
        lJcel=[]
        cells = {}
        
        #read first line first to get the first time point (could probably more clever, but hey)
        for line in fin:
            line=line.split(' ')
            prevtime=int(line[0])
            cells[int(line[1])]=int(line[11])
            break


        #read rest of file
        for line in fin:

            line=line.rstrip('\n').split(' ')
            time=int(line[0])
            thiscell=int(line[1])
            
            #################################################
            ## this is where you'll want to change things! ##
            #################################################
            
            if not ((prevtime+10000)%season): #almost at the end of the season: register adhesion values --> you might want to do this differently, or decide that "season" is just how often you sample adhesion values
                #print "reading jvals"
                if line[18]=='0': ## Check whether this is still the correct position!!!!! if it is '0' that means that this cell is in contact with the medium; the next position is the jvalue with the medium
                    lJmed.append(int(line[19])) #this is J val with med
                    if len(line)>20: # this cell may also be in contact with other cells, besides the medium; if you had to change the index above, also this has to change
                        lJcel.extend([ int(x) for x in line[21::2] ] ) 
                else:
                    lJcel.extend([ int(x) for x in line[19::2] ] ) 
            
            
            #were at a new time point
            if time!=prevtime:
                prevtime=time

		 #new season: put data in storage variables and reset vars
                if not(prevtime%season):
                    if (seasoncount>1 ): # or seasoncount<20:
                        gammlist=[med-cel/2. for med in lJmed for cel in lJcel]
                        if not gammlist:
                            gammlist = [0]
                        adval[0].append(np.median(gammlist))
                        storetime.append(seasoncount)  ### you might want to change this too!!!!

                    #reset vars
                    lJmed=[]
                    lJcel=[]
                    seasoncount+=1
                   
            

    #to plot gamma over time
    if (filecount==1): #set axis labels and limits
        ax0.set_ylabel("gamma",color="royalblue",fontsize=14)
        ax0.set_ylim(-6, 14)

        ax0.plot(storetime,adval[0],c=gamcolorset[filecount-1], lw=1.5, alpha=0.75)

plt.show()
plt.tight_layout()

fig.savefig(figname)
#plt.show()
