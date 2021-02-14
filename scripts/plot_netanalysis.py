import sys
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) <3:
    print "This is the program 'plot_netanalysis.py'"
    print "Usage: ./plot_division_time.py <input filename> <figure name> "
    sys.exit(1)
else:
    infile=sys.argv[1]
    outfile=sys.argv[2]
    minintercept=int(sys.argv[3])


divmig = np.loadtxt(infile, dtype='i', delimiter='\t')
sizes = np.shape(divmig[1:,1:])    
print sizes 

#fig, ax = plt.subplots(1,1)
fig=plt.figure() #!
fig.set_size_inches(1, 1.*sizes[0]/sizes[1], forward = False) #!
ax = plt.Axes(fig, [0., 0., 1., 1.]) #!
ax.set_axis_off() #!
fig.add_axes(ax) #!
ax.imshow(divmig[1:,1:], cmap='RdYlBu', origin='lower')

#find 2 parameters
pos=[0]*sizes[1]
found=[0]*sizes[1]
slope=0.
for i in range(sizes[0]):
    for j in range(sizes[1]):
        if divmig[1+i,j+1] == 1 and found[j]==0:
            pos[j]= i  #divmig[1+i,0]
             
            #print "pos ", pos[j]
            found[j]=1


#find the intercept and slope (for plotting you need the coordinates, not the actual gradient values)
if found[0]==0:
  intercept=divmig[-1,0]+minintercept
  plotintercept=sizes[0]
  slope=0
  plotslope=0
else:
    intercept=divmig[1+pos[0],0]
    plotintercept=pos[0]
    slope=0.
    plotslope=0.
    for ind, x in enumerate(pos):
        if ind >0 and found[ind] and found[ind-1]:
            slope+=divmig[1+x,0]-divmig[1+pos[ind-1],0]
            plotslope+=x-pos[ind-1]
        elif not found[ind]:
            if found[ind-1]:
                plotslope+=sizes[0]-pos[ind-1]
                slope+=divmig[-1,0]+minintercept-divmig[1+pos[ind-1],0]
            else:
             plotslope+=0
             slope+=0
    slope/=sizes[1]-1
    plotslope/=sizes[1]-1
    #ax.set_ylim(0,divmig[-1,0])
    #ax.plot(divmig[0,1:], plotintercept+plotslope*divmig[0,1:], color='black')

print "intercept ", intercept, " slope ", slope

#ax.set_xticks([ x for x in range(len(divmig[0,1:])) ] )
#ax.set_yticks([ x for x in range(len(divmig[1:,0])) ] )
#ax.set_xticklabels(divmig[0,1:])
#ax.set_yticklabels(divmig[1:,0])
#ax.set_xticklabels([])
#ax.set_yticklabels([])
#ax.set_xticks([])
#ax.set_yticks([])

divshare=divmig[1:,1:].sum()
migshare=(sizes[0])*(sizes[1])-divshare
writefile=open(outfile+"_pars.dat", "w")
writefile.write(str(intercept)+"\t"+str(slope)+"\n")
writefile.close()
migs="%04d" % (migshare,)
#print divs
plt.savefig("div_"+str(migs)+"_"+outfile+".pdf", dpi=sizes[1]) #bbox_inches='tight'
plt.close()

