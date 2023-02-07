#!/usr/bin/python2.7

'''
Version 2.
Plots pCellEvol data for the evolving moving cells model
'''

'''
minimal plotting tool for pCellEvol pred prey
the file is organised like this:
for all pop, every so many time steps (typically 1000), each individual is dumped, information is:
Time tau birthday key lock tau_contact J_contact tau_contact J_contact tau_contact J_contact ...
'''

import sys,math,os,subprocess
#from PIL import Image
import matplotlib as mpl
#mpl.use('Agg')
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
      maxwidht=50000
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
      color = 'orangered'
      #ax.plot(list_pos,lextra, color=facecolor,linestyle='-',marker='D',markersize=5, lw=1.5)
      ax.plot(list_pos,lextra, color=facecolor,linestyle='None',marker='D',markersize=5, lw=1.5)
      
#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################


filename=""
switchfile=""

maxTime=1000000000000000
if len(sys.argv)>1:
  pos=1
  while pos < len(sys.argv) :
    if sys.argv[pos]=='-maxtime':
      pos+=1
      maxTime=int(sys.argv[pos])
    elif sys.argv[pos]=='-filename':
      #print "Changing filename to",
      pos+=1
      filename=sys.argv[pos]
      print filename  
    elif sys.argv[pos]=='-switchfile':
      pos+=1
      switchfile=sys.argv[pos]
      print switchfile
    else:
      print "This is the program 'plot_data_cellcount.py'"
      print "It plots the birthday of pred and preys alive at a time point"
      print "Here are options, filename is necessary:"
      print "-filename [filename]"
      print "-maxtime [INT]"
      print "-switchfile [filename]"
      sys.exit(1)
    pos +=1

ldata=[]
ltime=[]
l_avrgdata_time=[]
lcounter=np.zeros((3,3),dtype=int)
lJmatrix=np.zeros((3,3),dtype=int)

l_avrgdata_time2=[]
lcounter2=np.zeros((3,3),dtype=int) #version 2 keeps count of each single value,so that proportions can be made for blob or scatter plot 
                                    # *** actually this might be useless, as blobplot funtion makes internal histogram

lJmatrix2 = [[[] for j in xrange(3)] for i in xrange(3)]

l_avrgdata_time_blob=[]
l_avrgdata_time2_blob=[]

lcounter_blob=np.zeros((3,3),dtype=int)
lJmatrix_blob=np.zeros((3,3),dtype=int)
lJmatrix2_blob = [[[] for j in xrange(3)] for i in xrange(3)]

#print lJmatrix2
#sys.exit(1)


lgamma_time=[]
lgamma=3*[0]  # this is G01, G02, G12
lgammapop=3*[0]
lpop=[0,0,0]
lpop_time=[]
lmaintf=[]
lmaintf_thistime=[]
lJfract=[]
lJfract_thistime=[]
lcfract=[]
lcfract_thistime=[]
lk=[]
lj=[]
lc=[]
lk_thistime=[]
lj_thistime=[]
lc_thistime=[]

#open switch time file and save time of switch:
switchtime=[]
if switchfile:
    with open(switchfile,"r") as fin:
      for line in fin:
        try:
          line=line.split()
          switchtime.append(int(line[-1]))
        except:
          pass

#open filename once to see what is the first time step... for now it is typically zero but in the future?
with open(filename,"r") as fin:
  for line in fin:
    line=line.split()
    inittime=int(line[0])
    break
ltime.append(inittime)

ldata=[]
# print "Hello0"
with open(filename,"r") as fin:
      # print "Hello1"
      for line in fin:
        # print line
        line=line.split()
        
        time=int(line[0])
        if time>maxTime: break
        sigma=int(line[1])
        tau=int(line[2])
        
        meanx = float(line[3])
        meany = float(line[4])
        vecx = float(line[5])
        vecy = float(line[6])
        chemx = float(line[7])
        chemy = float(line[8])
        
        birthdate=int(line[9])
                
        #key=line[10]  #don't care about this right now
        #lock=line[11]
        mu = float(line[12])
        particle = float(line[13])
        
        #print contacts
        #maintf=float(line[5])*(1-float(line[7])) #actually prints maintf
        
        maintf=float(line[14])  #fraction of resources to maintenance
        Jfract=float(line[19]) #fraction of own j val expressed
        cfract=[float(line[24])] #fraction of movement dedicated to gradient following
                
        km0,kmA,kmP,kmC = [ float(x) for x in line[15:19] ]  # maint vs movement
        kj0,kjA,kjP,kjC = [ float(x) for x in line[20:24] ] # expression of J val
        kc0,kcA,kcP,kcC = [ float(x) for x in line[25:29] ] # chemotaxis
        
        contacts=line[29:]
        # print contacts
        if time != ltime[-1]: 
          #calculate statistics
          lJ_av_matrix = np.divide(lJmatrix.astype(float) , lcounter.astype(float), out=np.zeros_like(lJmatrix.astype(float)), where=lcounter!=0)
          lJ_av_matrix_blob = np.divide(lJmatrix_blob.astype(float) , lcounter_blob.astype(float), out=np.zeros_like(lJmatrix_blob.astype(float)), where=lcounter_blob!=0)
          lk.append(lk_thistime) 
          lj.append(lj_thistime)          
          lc.append(lc_thistime)
          
          #lJmatrix.astype(float)/lcounter.astype(float)
          l_avrgdata_time.append(lJ_av_matrix)
          l_avrgdata_time2.append(lJmatrix2)
          
          l_avrgdata_time_blob.append(lJ_av_matrix_blob)
          l_avrgdata_time2_blob.append(lJmatrix2_blob)
          
          
          #append pop sizes
          lpop_time.append(lpop)
          #append the new time
          ltime.append(time)
          lmaintf.append( lmaintf_thistime )
          lJfract.append(lJfract_thistime)
          lcfract.append(lcfract_thistime)
          # zero the counters
          lcounter=np.zeros((3,3),dtype=int)
          lJmatrix=np.zeros((3,3),dtype=int)
          lJmatrix2=[[[] for j in xrange(3)] for i in xrange(3)]
          
          lcounter_blob=np.zeros((3,3),dtype=int)
          lJmatrix_blob=np.zeros((3,3),dtype=int)
          lJmatrix2_blob=[[[] for j in xrange(3)] for i in xrange(3)]
          
          lpop=[0 for _ in lpop]
          lmaintf_thistime=[]
          lJfract_thistime=[]
          lcfract_thistime=[]
          lk_thistime=[]
          lj_thistime=[]
          lc_thistime=[]
          
        lmaintf_thistime.append(maintf)
        lk_thistime.append( [ km0,kmA,kmP,kmC ] )
        lJfract_thistime.append(Jfract)
        lj_thistime.append( [ kj0,kjA,kjP,kjC ] )
        lcfract_thistime.append(cfract)
        lc_thistime.append([kc0,kcA,kcP,kcC])        
        
        lpop[tau]+=1 
        for ctau, cJ in [ contacts[pos:pos + 2] for pos in xrange(0, len(contacts), 2) ]:
          ctau=int(ctau)
          # print tau,ctau
          # print line
          # print contacts
          lcounter[tau,ctau]+=1 
          lJmatrix[tau,ctau]+=int(cJ)
          lJmatrix2[tau][ctau].append(int(cJ))
          if int(cJ)<=0: 
            print "Got cJ<=0 in line:", line
        
        # repeat above, but get only cells that are in contacts with more than just medium,
        # i.e. cells in a blob    
        if len(contacts) > 6:
          for ctau, cJ in [ contacts[pos:pos + 2] for pos in xrange(0, len(contacts), 2) ]:
            ctau=int(ctau)
            # print tau,ctau
            # print line
            # print contacts
            lcounter_blob[tau,ctau]+=1 
            lJmatrix_blob[tau,ctau]+=int(cJ)
            lJmatrix2_blob[tau][ctau].append(int(cJ))
            if int(cJ)<=0: 
              print "Got cJ<=0 in line:", line
        
        
print("Done reading file, now processing")
lpop_time=zip(*lpop_time)
ltime=ltime[:-1]  #removing last data point because it is not appended in the loop above

#print ltime

#print len(ltime)
#print len(lpop_time[0])
#print len(l_avrgdata_time)

#print [ J_av[1,2] for J_av in l_avrgdata_time ]

fig = plt.figure()
# ax1 = fig.add_subplot(411)
# ax1 = plt.subplot2grid((6, 2), (0, 1), colspan=3) #colspan makes it larger than a single column
ax1 = plt.subplot2grid((6, 4), (0, 0), colspan=4)
# print ltime
# print lpop_time[1]
# print "Hello"
ax1.plot( ltime,lpop_time[1], label="#prey" )
try: 
  ax1.plot( ltime,lpop_time[2], label="#pred" )
except:
  pass
if switchtime:
  ax1.scatter( switchtime, [1]*len(switchtime),marker="|" , label="switch" )
ax1.legend()
print("Hello1")

# ax2 = fig.add_subplot(412)
ax2 = plt.subplot2grid((6, 4), (1, 0), colspan=4)
#THIS BELOW WORKS FINE, BUT NOW WE USE GAMMAs, NOTICE THAT I AM USIGN AVERAGE J VALUES
a=ax2.plot( ltime, [ J_av[1,0] for J_av in l_avrgdata_time ], label="<J10>", lw = 0.5 )
b=ax2.plot( ltime, [ J_av[1,1] for J_av in l_avrgdata_time ], label="<J11>", lw = 0.5 )
c=ax2.plot( ltime, [ J_av[1,0] - 0.5*J_av[1,1] for J_av in l_avrgdata_time ], label="$\langle \langle \gamma \\rangle \\rangle$" )
#same for cells in blob
e=ax2.plot( ltime, [ J_av[1,0] for J_av in l_avrgdata_time_blob ], label="<J10>b", lw = 0.5 )
f=ax2.plot( ltime, [ J_av[1,1] for J_av in l_avrgdata_time_blob ], label="<J11>b", lw = 0.5 )
g=ax2.plot( ltime, [ J_av[1,0] - 0.5*J_av[1,1] for J_av in l_avrgdata_time_blob ], label="$\langle \langle \gamma b\\rangle \\rangle$" , lw = 0.5)

d=ax2.plot( ltime, [ 0 for _ in ltime ], lw = 0.5, linestyle='dashed')
# c=ax2.plot( ltime, [ J_av[1,2] for J_av in l_avrgdata_time ], label="J12" )
print("Av. last 10 gen. J10,J11 =", np.mean([ J_av[1,0] for J_av in l_avrgdata_time[-10:] ]), np.mean([ J_av[1,1] for J_av in l_avrgdata_time[-10:] ]) )
print("Av. last 10 gen. J10,J11 =", np.mean([ J_av[1,0] for J_av in l_avrgdata_time_blob[-10:] ]), np.mean([ J_av[1,1] for J_av in l_avrgdata_time_blob[-10:] ]) )

print("Av. last 10 gen. gamma", np.mean([ J_av[1,0] - 0.5*J_av[1,1] for J_av in l_avrgdata_time[-10:] ]))
print("Av. last 10 gen. gamma in blob", np.mean([ J_av[1,0] - 0.5*J_av[1,1] for J_av in l_avrgdata_time_blob[-10:] ]))
# d=ax2.plot( ltime, [ J_av[2,0] for J_av in l_avrgdata_time ], label="J20" )
#e=ax2.plot( ltime, [ J_av[2,1] for J_av in l_avrgdata_time ], label="J21" )  #overlaps perfectly with [1,2], as it should
# f=ax2.plot( ltime, [ J_av[2,2] for J_av in l_avrgdata_time ], label="J22" )

print("Hello1.1")
# print()
blobplot_every= 20
#blob_plot(ax2,[x[1][0] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (a[0].get_color(),0.5),'int', extra='median')
#blob_plot(ax2,[x[1][1] for x in l_avrgdata_time2][blobplot_every/2::blobplot_every] ,ltime[blobplot_every/2::blobplot_every], (b[0].get_color(),0.5),'int', extra='median')


#blob_plot(ax2,[x[1][2] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (c[0].get_color(),0.5),'int',extra='mean')
#blob_plot(ax2,[x[2][0] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (d[0].get_color(),0.5),'int')
#blob_plot(ax2,[x[2][1] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (e[0].get_color(),0.5),'int') #overlaps perfectly with [1,2], as it should
#blob_plot(ax2,[x[2][2] for x in l_avrgdata_time2][::blobplot_every] ,ltime[::blobplot_every], (f[0].get_color(),0.5),'int')
#print [x[2][2] for x in l_avrgdata_time2]

#blob_plot(ax2,[x[2][2] for x in l_avrgdata_time2][::10] ,ltime[::10], ('yellow',0.8),'int')

#plt.show()
#print ltime
#sys.exit(1)

#PLOT GAMMAS  
# ~G10~ = <J10> - <J11>/2
# ~G20~ = J20 - J22/2
# ~G12~ = J12 - (J11+J22)/2
          
#ax2.plot( ltime, [ J_av[1,0]-J_av[1,1]/2. for J_av in l_avrgdata_time ], label="G10" )
#ax2.plot( ltime, [ J_av[2,0]-J_av[2,2]/2. for J_av in l_avrgdata_time ], label="G20" )
#ax2.plot( ltime, [ J_av[1,2]- (J_av[1,1]+J_av[1,1])/2. for J_av in l_avrgdata_time ], label="G12" )


ax2.legend()
print("Hello2 - for this version Im not printing maint fract, only contact fract and chemotaxis fraction")
# ax3 = fig.add_subplot(413)
#sys.exit(1) #I have not impleneted it yet

# ax3 = plt.subplot2grid((6, 4), (2, 0), colspan=2)
# a = ax3.plot(ltime,[np.mean(x) for x in lmaintf],label="maint fract")
# blob_plot(ax3,lmaintf[::blobplot_every] ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
# ax3.legend()
# print("Hello3")
# 
# ax4 = plt.subplot2grid((6, 4), (3, 0))
# # at every time step lj is like this [[k0,kA,kP,kC],[k0,kA,kP,kC],...]
# a = ax4.plot( ltime, [ np.mean( zip(*x)[0] ) for x in lj ], label="km0")
# blob_plot(ax4, zip(*lj[::blobplot_every])[0] ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
# ax4.legend()
# 
# ax5 = plt.subplot2grid((6, 4), (3, 1))
# # at every time step lj is like this [[k0,kA,kP,kC],[k0,kA,kP,kC],...]
# a = ax5.plot( ltime, [ np.mean( zip(*x)[1] ) for x in lj ], label="kmA")
# blob_plot(ax5, zip(*lj[::blobplot_every])[1] ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
# ax5.legend()
# 
# ax6 = plt.subplot2grid((6, 4), (4, 0))
# # at every time step lj is like this [[k0,kA,kP,kC],[k0,kA,kP,kC],...]
# a = ax6.plot( ltime, [ np.mean( zip(*x)[2] ) for x in lj ], label="kmP")
# blob_plot(ax6, zip(*lj[::blobplot_every])[2] ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
# ax6.legend()
# 
# ax7 = plt.subplot2grid((6, 4), (4, 1))
# # at every time step lj is like this [[k0,kA,kP,kC],[k0,kA,kP,kC],...]
# a = ax7.plot( ltime, [ np.mean( zip(*x)[3] ) for x in lj ], label="kmC")
# blob_plot(ax7, zip(*lj[::blobplot_every])[3] ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
# ax7.legend()

ax3 = plt.subplot2grid((6, 4), (2, 0), colspan=2)
a = ax3.plot(ltime,[np.mean(x) for x in lcfract],label="chemot fract")
blob_plot(ax3,lcfract[::blobplot_every] ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
ax3.legend()
print("Hello3")

lc_toblpl=lc[::blobplot_every] # same as before, just sparser

ax4 = plt.subplot2grid((6, 4), (3, 0))
# at every time step lj is like this [[k0,kA,kP,kC],[k0,kA,kP,kC],...]
a = ax4.plot( ltime, [ np.mean( zip(*x)[0] ) for x in lc ], label="kc0")


bla = []
for thistime_lc_toblpl in lc_toblpl:
    transp=zip(*thistime_lc_toblpl)
    bla.append(transp[0])
blob_plot(ax4, bla ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')    
#blob_plot(ax4, zip(*lc[::blobplot_every])[0] ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
ax4.legend()

ax5 = plt.subplot2grid((6, 4), (3, 1))
# at every time step lj is like this [[k0,kA,kP,kC],[k0,kA,kP,kC],...]
a = ax5.plot( ltime, [ np.mean( zip(*x)[1] ) for x in lc ], label="kcA")
bla = []
for thistime_lc_toblpl in lc_toblpl:
    transp=zip(*thistime_lc_toblpl)
    bla.append(transp[1])
blob_plot(ax5, bla ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
ax5.legend()

ax6 = plt.subplot2grid((6, 4), (4, 0))
# at every time step lj is like this [[k0,kA,kP,kC],[k0,kA,kP,kC],...]
a = ax6.plot( ltime, [ np.mean( zip(*x)[2] ) for x in lc ], label="kcP")
bla = []
for thistime_lc_toblpl in lc_toblpl:
    transp=zip(*thistime_lc_toblpl)
    bla.append(transp[2])
blob_plot(ax6, bla ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
ax6.legend()

ax7 = plt.subplot2grid((6, 4), (4, 1))
# at every time step lj is like this [[k0,kA,kP,kC],[k0,kA,kP,kC],...]
a = ax7.plot( ltime, [ np.mean( zip(*x)[3] ) for x in lc ], label="kcC")
bla = []
for thistime_lc_toblpl in lc_toblpl:
    transp=zip(*thistime_lc_toblpl)
    bla.append(transp[3])
blob_plot(ax7, bla ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
ax7.legend()


ax8 = plt.subplot2grid((6, 4), (2, 2), colspan=2)
# ax4 = fig.add_subplot(414)
a = ax8.plot(ltime,[np.mean(x) for x in lJfract],label="J expr fract")
blob_plot(ax8,lJfract[::blobplot_every] ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
ax8.legend()

lj_toblpl=lj[::blobplot_every] # same as before, just sparser

ax9 = plt.subplot2grid((6, 4), (3, 2))
# at every time step lj is like this [[k0,kA,kP,kC],[k0,kA,kP,kC],...]
a = ax9.plot( ltime, [ np.mean( zip(*x)[0] ) for x in lj ], label="kj0")
bla = []
for thistime_lj_toblpl in lj_toblpl:
    transp=zip(*thistime_lj_toblpl)
    bla.append(transp[0])
blob_plot(ax9, bla ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
ax9.legend()

ax10 = plt.subplot2grid((6, 4), (3, 3))
# at every time step lj is like this [[k0,kA,kP,kC],[k0,kA,kP,kC],...]
a = ax10.plot( ltime, [ np.mean( zip(*x)[1] ) for x in lj ], label="kjA")
bla = []
for thistime_lj_toblpl in lj_toblpl:
    transp=zip(*thistime_lj_toblpl)
    bla.append(transp[1])
blob_plot(ax10, bla ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
ax10.legend()

ax11 = plt.subplot2grid((6, 4), (4, 2))
# at every time step lj is like this [[k0,kA,kP,kC],[k0,kA,kP,kC],...]
a = ax11.plot( ltime, [ np.mean( zip(*x)[2] ) for x in lj ], label="kjP")
bla = []
for thistime_lj_toblpl in lj_toblpl:
    transp=zip(*thistime_lj_toblpl)
    bla.append(transp[2])
blob_plot(ax11, bla ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
ax11.legend()

ax12 = plt.subplot2grid((6, 4), (4, 3))
# at every time step lj is like this [[k0,kA,kP,kC],[k0,kA,kP,kC],...]
a = ax12.plot( ltime, [ np.mean( zip(*x)[3] ) for x in lj ], label="kjC")
bla = []
for thistime_lj_toblpl in lj_toblpl:
    transp=zip(*thistime_lj_toblpl)
    bla.append(transp[3])
blob_plot(ax12, bla ,ltime[::blobplot_every], (a[0].get_color(),0.5),'float')
ax12.legend()

title = "no title"
try:
  title=filename.split('/')[-3]+'_'+filename.split('/')[-2]+'_'+filename.split('/')[-1]
except:
  pass
  
fig.suptitle(title)

# Option 2
# TkAgg backend
manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
#plt.savefig(title+'.pdf')
#ax2.set_xlim([0,50000])

plt.savefig("bla.pdf",bbox_inches='tight')
plt.show()
sys.exit(1)
