#!/usr/bin/python2.7

'''
Plot the flow field, which is a fancy version of the equal time correlation resolved in space
To use only on one file at a time, possibly with larga gamma, if you want decent results.

The idea is as follows:
for every cell, calculate its displacement as a vector, 
then consider every cell as the center and for every other cell also calculate its dipl. relative to the first cell, that is a fter rotating second cell by th angle of the first cell.
At this point, bin in space and take the average (rotated) vector.
 
'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
# mpl.use('Agg')
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
# fig, (ax0,ax1) = plt.subplots(nrows=2,sharex=True)
# fig.subplots_adjust(hspace=0)

#fig, ax0 = plt.subplots()

if len(sys.argv) <3:
  print "This is the program 'plot_flowfield_fromCM.py'"
  print "Usage: ./plot_flowfield_fromCM.py <figure name> <cell targetarea> <binsize [pix]> <filename> <filename> <filename> ..."
  sys.exit(1)
else:
  figname=sys.argv[1]
  cell_tarsize=float(sys.argv[2])
  binsize=int(sys.argv[3])
  
macoldist=0
filecounter=0


field_size = 2*500

binsize=10 #8 is roughly 2 radiuses
dimensions_flowfield = field_size/binsize 
# flow_field = np.array((512,512), dtype=[0,0,0])
flow_col = np.zeros((dimensions_flowfield,dimensions_flowfield), dtype= float)
flow_row = np.zeros((dimensions_flowfield,dimensions_flowfield), dtype= float)
flow_howmany = np.zeros((dimensions_flowfield,dimensions_flowfield), dtype= float)



for filename in sys.argv[4:]:
  colpos=[]
  rowpos=[]
  xn=[]
  colpos.append([])
  rowpos.append([])
  count=0

  print filename
  # avangle=[0.]*200*int(binmult)   #1 bin is binsize celldiam 
  # binamount=[0]*200*int(binmult)  #to keep track of how many in each bin
  # angles=[[] for x in range(200*int(binmult))]
  
  ##################################################
  ##read data from file: store raw cell positions##
  ##################################################
  
  #read time in the first line
  output = subprocess.Popen(['head', '-1', filename], stdout=subprocess.PIPE).communicate()[0]
  time0 = int(output.split(' ')[0])
  timepoint=time0
  #read file
  with open(filename,"r") as fin:
    for line in fin:
      line=line.split()
      if int(line[0])!=timepoint:
        colpos.append([-1.]*len(colpos[0]))
        rowpos.append([-1.]*len(rowpos[0]))
        count=0
        timepoint=int(line[0])
      if timepoint==time0:
          # if int(line[31])==0:
          #   rowpos[-1].append(-1.)
          #   colpos[-1].append(-1.)
          # else:
            rowpos[-1].append(float(line[3]))
            colpos[-1].append(float(line[4]))
      else:
          # if int(line[31])!=0:
              rowpos[-1][int(line[1])-1]=float(line[3])
              colpos[-1][int(line[1])-1]=float(line[4])
          # else:
          #     rowpos[-1][int(line[1])-1]=-1
          #     colpos[-1][int(line[1])-1]=-1
      count+=1
  
  # Center of mass
  row_cm = [np.mean(rows) for rows in rowpos]
  col_cm = [np.mean(cols) for cols in colpos]
  
  maxint=len(colpos)
  nrcells=len(colpos[0])
  print "maxint=",maxint
  deltapos = 5 # when we do x_later - x_now, delta says how many positions back we look
  for i in range(500,maxint,100): #skip first time steps, not to include initial conditions and also so that i>deltapos
    if i*50<deltapos: continue
    count=0
    count2=0
    xn=[]      
    
    flow_center_row =  int(0.5+(field_size/2)/float(binsize))
    flow_center_col = flow_center_row
    
    row_cm_mov = row_cm[i] - row_cm[i-deltapos] #displacement of the center of mass
    col_cm_mov = col_cm[i] - col_cm[i-deltapos]
    modulus_cm_mov1 = math.sqrt(row_cm_mov**2 + col_cm_mov**2)
    if modulus_cm_mov1<0.0001: continue
    cos_phi = col_cm_mov / modulus_cm_mov1 # for cells going left this should be about -1
    sin_phi = row_cm_mov / modulus_cm_mov1 # and this 0
    
    # alternatively I could use as a preferred frame of reference 
    # the direction from center of mass to peak of the gradient
    
    # flow_row[ flow_center_row , flow_center_col ] += row_cm_mov
    # flow_col[ flow_center_row , flow_center_col ] += col_cm_mov
    # flow_howmany[flow_center_row,flow_center_col ] += 1
    # 
    for c2 in range(nrcells):
        if colpos[i][c2]>-1 and colpos[i-deltapos][c2]>-1:
            #Delta movement and length of dsiplacement
            
            #dist=math.sqrt((colpos[i][c1]-colpos[i][c2])**2+(rowpos[i][c1]-rowpos[i][c2])**2) # this is the distance after the movement
            colmov2=colpos[i][c2]-colpos[i-deltapos][c2] #displacement of cell c2
            rowmov2=rowpos[i][c2]-rowpos[i-deltapos][c2]

            #rotate cell 2 by the angle of cell one, this is done by the rotation matrix
            # ( cos(phi) -sin(phi) )
            # ( sin(phi)  cos(phi) )
            # now, cos(phi) = colmov1 / dist (length of the x component of vector scaled by total length)
            # and sin(phi) = rowmov1 / dist
            colmov2_rot =   colmov2 * cos_phi + rowmov2 * sin_phi
            #colmov2_rot = colmov2 * sin_phi - rowmov2 * cos_phi
            rowmov2_rot = - colmov2 * sin_phi + rowmov2 * cos_phi 

            # print "colmov1,rowmov1, phi =",colmov1,rowmov1, 360*math.atan2( cos_phi,sin_phi)/(2.*math.pi) 
            # print "colmov2,rowmov2 =",colmov2,rowmov2
            # print "colmov2_rot,rowmov2_rot =",colmov2_rot,rowmov2_rot
            # 


            #x-distance and y-distance between center of mass and c2 center of mass, at time i-deltapos
            rowdist = rowpos[i-deltapos][c2] - row_cm[i-deltapos] #+ field_size/2
            coldist = colpos[i-deltapos][c2] - col_cm[i-deltapos] #+ field_size/2 # to center it
            
            print "rowdist,coldist",rowdist,coldist
            
            coldist_rot = coldist * cos_phi + rowdist * sin_phi # + field_size/2
            rowdist_rot = - coldist * sin_phi + rowdist * cos_phi  #+ field_size/2

            dist = np.hypot(coldist,rowdist)  #distance c1,c2
            dist_rot = np.hypot(coldist_rot,rowdist_rot) #distace rotated c2 from rotated c1
            if dist > dist_rot+0.01 or dist < dist_rot- 0.01 or dist < 4 or dist_rot < 4: 
              print "Error:"
              print "Time: ", i*20
              # print "c1, c2: ", c1,c2
              print "dist, dist_rot", dist, dist_rot
              # print "rowpos 1, colpos 1", rowpos[i-deltapos][c1], colpos[i-deltapos][c1]
              # print "rowpos 2, colpos 2", rowpos[i-deltapos][c2], colpos[i-deltapos][c2]
              
            # else:
            #     print "Hello"
            flow_rowbin = int(0.5+(rowdist_rot+ field_size/2)/float(binsize))
            flow_colbin = int(0.5+(coldist_rot+ field_size/2)/float(binsize))
            flow_row[       flow_rowbin , flow_colbin ] += rowmov2_rot
            flow_col[       flow_rowbin , flow_colbin ] += colmov2_rot
            flow_howmany[ flow_rowbin , flow_colbin ] += 1
      # sys.exit(1)  
print "Done operations, now plotting"
print "x,y for printing :",(field_size/2)/binsize
print flow_col[(field_size/2)/binsize, (field_size/2)/binsize] 
print flow_row[(field_size/2)/binsize, (field_size/2)/binsize]
print flow_howmany[(field_size/2)/binsize, (field_size/2)/binsize]
# sys.exit(1)

# flow_col = np.divide(flow_col,flow_howmany, out=np.zeros_like(flow_col), where=flow_howmany!=0) #returns zero where it divides by zero, because there is no flow
# flow_row = np.divide(flow_row,flow_howmany, out=np.zeros_like(flow_row), where=flow_howmany!=0)
flow_col = np.divide(flow_col,flow_howmany, out=np.zeros_like(flow_col), where=flow_howmany>1*len(sys.argv[4:])) #returns zero where it divides by zero, because there is no flow
flow_row = np.divide(flow_row,flow_howmany, out=np.zeros_like(flow_row), where=flow_howmany>1*len(sys.argv[4:])) #also returns zero if there aren't enough values to make an evarge

# for r,c in np.ndindex(flow_col.shape):
#   if flow_col[r,c]<0.00001 and flow_row[r,c]<0.00001:
#       flow_col[r,c]=np.nan
#       flow_row[r,c]=np.nan
# 
r = np.linspace(0,dimensions_flowfield,dimensions_flowfield)
c = np.linspace(0,dimensions_flowfield,dimensions_flowfield)

C,R = np.meshgrid(c,r)

print "After some filtering: "
print flow_col[(field_size/2)/binsize, (field_size/2)/binsize] 
print flow_row[(field_size/2)/binsize, (field_size/2)/binsize]

fig, ax = plt.subplots()
if deltapos!=5:    
    print "WAAAARNING: Deltapos is not 5"

bla= ax.quiver(C,R, flow_col, flow_row,  flow_howmany, cmap=plt.get_cmap('viridis_r'), scale = 5*float(deltapos)/(binsize),angles='xy',scale_units='xy',  norm=mpl.colors.LogNorm())
fig.colorbar(bla)
# ax.quiver(X,Y, flow_col, flow_row, scale = 1)
# ax.set_xlim([(field_size/3)/binsize , (2*field_size/3)/binsize])
# ax.set_ylim([(field_size/3)/binsize , (2*field_size/3)/binsize])
ax.set_xlim([43 , 58])
ax.set_ylim([43 , 58])
ax.set_aspect('equal')
ax.set_title('deltapos ='+str(deltapos))
a = ax.get_xticks().tolist()
ax.set_xticklabels([int(x) -50 for x in a], ha="left")

plt.show()
sys.exit(1)
  
  
  
  