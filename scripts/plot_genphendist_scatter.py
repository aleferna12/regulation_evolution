#!/usr/bin/python3.8
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import os

colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]

#this script has to make a scatter plot of genotypic vs phenotypic distance of evolved individuals to original network

#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

fig, ax = plt.subplots()

if len(sys.argv) <2:
    print ("This is the program 'plot_genphendist_scatter.py'")
    print ("Usage: ./plot_genphendist_scatter.py <figure name> < input file name 1> ...>")
    sys.exit(1)
else:
    figname=sys.argv[1]
    colorset=sns.color_palette("twilight",n_colors=6)
    #print ("{}".format(colorset))

#sys.exit()
count=0

#print masked
for filename in sys.argv[2:]:
    #decide the colour for the dots (depends on type and nr of the simulation)

    gen, var, phen = np.loadtxt(filename,unpack=True, usecols=[0,1,2])
    sizes=10000*var
    pathname=filename.split("_")
    ax.scatter(gen, phen, s=sizes, c=[colorset[count]], label=pathname[1])
    count+=1

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlabel("genotypic distance")
ax.set_ylabel("phenotypic distance")
ax.set_aspect('auto')                                                              
fig.savefig(figname)
plt.show()


