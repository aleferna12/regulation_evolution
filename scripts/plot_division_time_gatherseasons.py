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

if len(sys.argv) <4:
    print ("This is the program 'plot_division_time_gatherseasons.py'")
    print ("Usage: ./plot_division_time.py <figure name> <season duration> <nr seasons> <label1 filename1> ...")
    sys.exit(1)
else:
    figname=sys.argv[1]
    season=int(sys.argv[2])
    nrseasons=int(sys.argv[3])
    llabel=sys.argv[4::2] #wow I did not know this! So handy....
    lfilename=sys.argv[5::2]
    #end=int(sys.argv[4])
    #filename=sys.argv[5]

pos=[1,2,3]
#print (llabel)
#print (lfilename)
##read data from file: store division data
for label,filename in zip(llabel,lfilename):
    #we get the original division time from label (or just an identifying number)
    #last saved time step is
   
    output = subprocess.Popen(['tail', '-1', filename], stdout=subprocess.PIPE).communicate()[0]
    arrout=output.decode('ascii').split(' ')
    lasttime = int(arrout[0])
    if(nrseasons>1):
        finalseason = lasttime - lasttime%season #finds last completed season (data will be empty at that point)
        end = finalseason - 10000 #at the round times, the new season has already started (so cells have been cleared) -- so reduce by one save step
        start = finalseason - nrseasons*season
    else:
        end=lasttime
        start=0

    #print (label, filename)
    with open(filename,"r") as fin:

        #read first line first to get the first time point (could probably more clever, but hey)
        divtime=[]
        divtime.extend([] for i in range(3))
        seasoncount=1
        far=0
        add=0
        cells=[0 for x in range(20000)]
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
            if(thiscell>20000):
                print("oh shit")
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
                    cells=[0 for x in range(20000)]
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
    print ("{} {}".format(label,sum(divtime[0])/float(len(divtime[0])) ) )
    div1=ax0.violinplot(divtime,positions=pos,widths=0.5, showmedians=True, showextrema=False)
    #div1['bodies'].set_facecolor(colours[:3])
    div1['cmedians'].set_color([ colors.to_rgba(x) for x in colours[:3]  ]  )
    for ind,pc in enumerate(div1['bodies']):
        pc.set_facecolor(colours[ind])
    
    pos=[x+4 for x in pos]

xtickpos=[2+4*i for i in range(len(llabel))]
ax0.set_xticklabels(llabel)
ax0.set_xticks(xtickpos)

ax0.set_ylim(0,season)
ax0.set_xlabel('simulation')
ax0.set_ylabel('division time')

fig.savefig(figname, bbox_inches='tight')
#plt.show()
