#!/usr/bin/python2.7

'''
Makes a set of blob plots with the last so many time steps for a series of files (ordered by season, maybe?)
Usage: 
./plot_steadystate..blabla.py <output file name> <#timesteps, from the last> file1 file2 file3 ...
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
### --- BLOB PLOT --- ###
###                   ###
#########################
def blob_plot(ax,data,list_pos, colour,datatype,extra=False):
    lextra=[]
    for d,p in zip(data,list_pos):
      dowecolor=True
      if len(d)==0: 
          continue
      if not d: 
        dowecolor=False
        continue
      #this makes the histogram that constitutes the blob plot
      m=min(d) # the largest value, upper edge of the histogram
      M=max(d) # the lowest value, lower...
      
      #print p
      #his=np.bincount(d)
      #print his
      
      if datatype=='int':
        his, bins = np.histogram(d, bins=[ 5*x for x in xrange(2+ M/5 ) ])
        #if len(d)>4:
          #print bins
          #print his
          #x = (bins[:-1] + bins[1:]) / 2.  #centers 
          #print x
          #sys.exit(1)
        
        #x = (bins[:-1] + bins[1:]) / 2.  #centers 
        x = bins[:-1] # seems to work better than centers
        
        #his=np.trim_zeros(his)
        #his=np.append(his,[0])
        #his=np.append([0],his)
        #x=np.linspace(m-1,M+1,M-m+3 )
        
      elif datatype=='float':
        nbins=43
        x = np.linspace(m,M,nbins) # support for histogram of floats
        
        his,bins = np.histogram(d, bins=nbins)
      
      #maxwidht=0.15
      maxwidht=2500
      #scale_blob=1.
      max_his=max(his)
      #print max_his
      if max_his>0.:
        scale_blob=maxwidht/float(max_his)
        #scale_blob=10
      else: 
        scale_blob=0.
        print "Warning, max_his is 0." 
        
      shift_his_plus_pos =  [ p + h*scale_blob  for h in his]
      shift_his_minus_pos = [ p - h*scale_blob  for h in his]
            
      color_alpha_blobs=colour
      facecolor,alpha=color_alpha_blobs # assign color and transparency
      #this is the matplotlib function that does the trick
      
      #print "len x:",len(x), "len shift_his_minus_pos", len(shift_his_minus_pos)
      ax.fill_betweenx(x,shift_his_minus_pos, shift_his_plus_pos, linewidth=0., facecolor= facecolor, alpha=alpha)
      #calculates median or mean, if you want
      if extra=='median': lextra.append( np.median(d) )
      elif extra=='mean': lextra.append( np.mean(d) )
    #and plots it
    if extra != False and dowecolor==True:
      color = 'darkgoldenrod'
      #ax.plot(list_pos,lextra, color=facecolor,linestyle='-',marker='D',markersize=5, lw=1.5)
      ax.plot(list_pos,lextra, color=facecolor,linestyle='None',marker='D',markersize=5, lw=1.5, zorder=5)
      return lextra
    

#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""
fig, ax = plt.subplots()
if len(sys.argv) < 3:
  print "This is the program 'plot_steadysate...blabla.py'"
  print "./plot_steadystate..blabla.py <output file name> <#timesteps from the last> <label1> <file1> <label2> <file2> <label3> <file3>  ..."
  sys.exit(1)
else:
  figname=sys.argv[1]
  timesteps = int(sys.argv[2])
  llabel=sys.argv[3::2]
  lfilename=sys.argv[4::2]

if len(llabel) != len(lfilename): 
  print "Error, different number of labels and filenames?"
  print llabel
  print lfilename
  sys.exit(1)

#the figure
fig, ax = plt.subplots()
llJmed=[]
llJcel=[]
lmeanJc=[]
lmeanJm=[]
print "llabels", llabel
print "filenames", lfilename

#some smart looping to get season and filename
for label,filename in zip(llabel,lfilename):
    #we get season from label
    season = int(label)
    #last saved time step is
    output = subprocess.Popen(['tail', '-1', filename], stdout=subprocess.PIPE).communicate()[0]
    lasttime = int(output.split(' ')[0])
    # by the way, data is saved every 10000
    # so if the season is 10000 or less we take the last 10 time steps and that's it
    
    # if season<=10000:
    #     ltime = [lasttime - i*10000 for i in range(10)]
    # 
    # #else, we should take a time point for each of the last 10 seasons
    # # which one? ... the middle one!
    # # BUT WHAT IF SEASON / 2 is not divisible by 10000? -> should get closest point to that...
    # # so instead of middle one we do trimmed_lasttime - 100000
    # else:
    trimmed_lasttime = lasttime - lasttime%season 
    timepoint = trimmed_lasttime - 10000
    ltime = [lasttime - i*season for i in range(10)]
    # print ltime
    # Now we go get the date
    lJmed=[]
    lJcel=[]
    with open(filename,"r") as fin:
        print "Opening", filename
        print "ltime", ltime
        # for line in fin[::-1]:
        for line in reversed(fin.readlines()):
            # print line
            line=line.split(' ')
            time = int(line[0])
            # print time
            if time in ltime:
                # print "Hello"
                # print line
                #save lJo, lJ1
                if line[29]=='0': 
                    lJmed.append(int(line[30])) #the line[30] is J val with med
                    lJcel.extend([ int(x) for x in line[32::2] ] ) 
                else:
                    lJcel.extend([ int(x) for x in line[30::2] ] ) 
    
    llJcel.append(lJcel)
    llJmed.append(lJmed)
    # At this point we can start filling up our figure
    # meanJc=blob_plot(ax,[lJcel],[season], ('royalblue',1.),'int' ,extra='mean')
    # meanJm=blob_plot(ax,[lJmed],[season], ('firebrick',1.),'int' ,extra='mean')
    
    # lmeanJc.extend(meanJc)
    # lmeanJm.extend(meanJm)
    # ax.boxplot(lJmed, positions=[math.log(season)])
    # ax.boxplot(lJcel, positions=[math.log(season)])

#sets flier properties    
#flierprops = dict(marker='.', markerfacecolor='k', markersize=1,linestyle='none') <- this works perfectly, I just don't want to show fliers any more, because nothing Changes

position_offset = 1000 # so that Jc and Jm don't overlap
labels=[int(label[:-3]) for label in llabel ]
print "labels:", labels

ax.boxplot(llJcel, positions=[int(x)-2.5*position_offset if i==llabel.index(x) else int(x)+1.25*position_offset  for i,x in enumerate(llabel)], widths=(2000), showfliers=False, medianprops=dict(color='royalblue') )
lmeanJc= [ np.mean(lJcel) for lJcel in llJcel]
# ax.plot(llabel,lmeanJc,color='royalblue')
ax.plot(llabel,lmeanJc,linestyle = 'None',marker='o', markerfacecolor='None', markeredgecolor='royalblue', markersize = 1)

ax.boxplot(llJmed, positions=[int(x)-1.25*position_offset if i==llabel.index(x) else int(x)+2.5*position_offset  for i,x in enumerate(llabel)], widths=(2000), showfliers=False, medianprops=dict(color='darkgoldenrod') )
lmeanJm= [ np.mean(lJm) for lJm in llJmed]
# ax.plot(llabel,lmeanJm,color='darkgoldenrod')
ax.plot(llabel,lmeanJm,linestyle = 'None',marker='o', markerfacecolor='None', markeredgecolor='darkgoldenrod', markersize = 1)

# ax.plot(llabel , [ jm-jc/2. for jc,jm in zip(lmeanJc,lmeanJm)], color='seagreen')
ax.plot(llabel , [ jm-jc/2. for jc,jm in zip(lmeanJc,lmeanJm)], linestyle = 'None',marker='o', markerfacecolor='None', markeredgecolor='seagreen', markersize = 1)
ax.plot(llabel,[0. for _ in llabel], lw = 0.5, color='seagreen')

ax.set_xticklabels(labels)
ax.set_xticks([int(label) for label in llabel ])

ax.set_xlabel('$\\tau_s$ [$\cdot 10^3$ mcs]', fontsize=12)

ax.set_xlim([0,160000])

# plt.show()    
plt.savefig(figname)

sys.exit(1)
        
"""




      #read first line first to get the first time point (could probably more clever, but hey)
      for line in fin:
        line=line.split()
        timepoint=int(line[0])
        startxpos+=float(line[3])
        startypos+=float(line[4])
        count+=1
        break

      #read rest of file, to find positions at beginning and end of period
      for line in fin:
        line=line.split()
        if int(line[0])==timepoint:
          startxpos+=float(line[3])
          startypos+=float(line[4])
          count+=1
        elif int(line[0])==timepoint+period:
          endxpos+=float(line[3])
          endypos+=float(line[4])
          count2+=1

    #calculate average positions at beginning and end of period
    startxpos/=float(count)
    startypos/=float(count)
    endxpos/=float(count2)
    endypos/=float(count2)

    xx=endxpos-startxpos
    yy=endypos-startypos
    
    displace.append(math.sqrt(xx*xx+yy*yy))
    
 
##start plotting##

ax.scatter(offset,displace, c=colourset,alpha=0.5)
ax.set_xlim([0,totalcounter+1])
#custom xlabels
my_labels=ax.get_xticks().tolist()
counter=0
newlab=[]
for el in my_labels:
  if el.is_integer() and counter<len(labelset) and el>0.001:
    newlab.append(labelset[counter])
    counter+=1
  else:
    newlab.append('')

ax.set_xticklabels(newlab)
#print my_labels
#print labelset
#print newlab

plt.xlabel('simulation')
ylab='displacement in time period '+str(period)
plt.ylabel (ylab)
fig.savefig(figname, bbox_inches='tight')
#plt.show()




"""
