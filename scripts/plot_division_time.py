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

def violin_plot(ax,data1,pos, colour):
    '''
    create violin plots on an axis
   
    #dist = max(pos)-min(pos)
    #w = min(0.15*max(dist,1.0),0.5)
    lmedian=[]
    lmean=[]
    ''' 
    m1=int(min(data1))
    M1=int(max(data1))
    #print "violin: ", m1, M1
    nbins1=(M1-m1+1)/10000
    
    if not nbins1:
        nbins1=1
    x1=np.linspace(m1,M1,nbins1)
    his1,bins1 = np.histogram(data1, bins=nbins1)
    '''
    m2=int(min(data2))
    M2=int(max(data2))
    nbins2=M2-m2+1
    x2=np.linspace(m2,M2,nbins2)
    his2,bins2 = np.histogram(data2, bins=nbins2)
    '''
    #x1=np.linspace(0.,1.0,nbins1) #-1./(2.*nbins1)
    #his1,bins1 = np.histogram(data1, bins=nbins1, range=(0.,1.0))
    data1.sort()
    med=data1[len(data1)/2]      
    
    #print x1
    #print his1
    #shift and scale his to assign proper pos and size to fit in the space
    scale_violin=0.8*every
    shift_his_minus_pos = [pos- scale_violin*h/float(len(data1)) for h in his1]
    shift_his_plus_pos = [pos + scale_violin*h/float(len(data1)) for h in his1]
    
    #print len(his1)
    #print len(x1)
    
    ax.fill_betweenx(x1,shift_his_minus_pos, shift_his_plus_pos, linewidth=0.0, facecolor=colours[colour], alpha=0.5)
    ax.scatter(pos,med,s=8, c=colours[colour])  
  


colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
#fig, (ax0, ax1) = plt.subplots(nrows=2)
fig, ax0 = plt.subplots()

if len(sys.argv) <3:
    print "This is the program 'plot_division_time.py'"
    print "Usage: ./plot_division_time.py <figure name> <season duration> <plotinterval> filename"
    sys.exit(1)
else:
    figname=sys.argv[1]
    season=int(sys.argv[2])
    every=int(sys.argv[3])
    filename=sys.argv[4]

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
      cells[int(line[1])]=int(line[11])
      break

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
                    print "oh oh...", divnum, prevtime
                    sys.exit(1)
                #do some plotting
                if not (seasoncount%every): # or seasoncount<20:
                    print "plotting season ", seasoncount
                    if seasoncount<20:
                        addpos=prevtime/season*every
                    else:
                        addpos=prevtime/season
                    div1=ax0.violinplot(divtime,positions=[prevtime/season-0.2*every,prevtime/season, prevtime/season+0.2*every],widths=0.5*every, showmedians=True, showextrema=False)
                    #div1['bodies'].set_facecolor(colours[:3])
                    div1['cmedians'].set_color([ colors.to_rgba(x) for x in colours[:3]  ]  )
		    for ind,pc in enumerate(div1['bodies']):
			pc.set_facecolor(colours[ind])
                        #med=divtime[ind][len(divtime[ind])/2]
                        #ax0.scatter(20+addpos-0.2*every+0.2*every*ind,med,s=8, c=colours[ind])
                    #violin_plot(ax0,divtime[0][1:],20+prevtime/season-0.2*every, 0)
                    #violin_plot(ax0,divtime[1][1:],20+prevtime/season, 1)
                    #violin_plot(ax0,divtime[2][1:],20+prevtime/season+0.2*every, 2)
                    if(far):
                        ax0.scatter(prevtime/season, 0, c="black", s=8)
            
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
           



ax0.set_xlabel('season')
ax0.set_ylabel('division time')
ax0.set_title('timing of cell divisions in a season')

fig.savefig(figname, bbox_inches='tight')
#plt.show()
