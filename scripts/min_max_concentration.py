#!/usr/bin/python2.7

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

if len(sys.argv) <1:
    print "This is the program 'min_max_concentration.py'"
    print "Usage: ./min_max_concentration.py filename"
    sys.exit(1)
else:
    filename=sys.argv[1]


##read data from file: store division data
with open(filename,"r") as fin:


    minconc=999999
    maxconc=0
    lncount=0
    for line in fin:
      line=line.split(' ')
      lncount+=1
      if int(line[12]) < 0 or int(line[12]) >10000:
        print "ai ai ", int(line[12]), " cell ", line[1], " type ", line[2]
      else:
        if int(line[12])<minconc:
          minconc=int(line[12])
        if int(line[12])>maxconc:
          maxconc=int(line[12])
 
print minconc, maxconc
