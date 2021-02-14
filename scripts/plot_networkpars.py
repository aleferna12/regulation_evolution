
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from scipy import stats
from matplotlib.lines import Line2D
import numpy as np
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

#colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]

colours=["dimgray", "red"]
MINGRAD=50
MAXGRAD=1550

def density_estimation(m1, m2,min, max):
    xmin=min[0]-50
    xmax=max[0]
    ymin=min[1]
    ymax=max[1]
    X, Y = np.mgrid[xmin:xmax:80j, ymin:ymax:80j]                                                     
    #print X
    #print Y
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])                                                                        
    kernel = stats.gaussian_kde(values) #,bw_method=0.1                                                                
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z


def makemask_old(minx, maxx, miny, maxy):
    #calculate max slopes for mask
    maxmax=int((MAXGRAD-MINGRAD)/3.)
    slopesteps=60
    masking=np.zeros( (slopesteps+1, int((maxx)/20)+1)) #-MINGRAD
    #print masking
    #print maxmax
    #print "looping..."
    for ind1, grad in enumerate(range(int(minx), int(maxx), 20)):
        #print "grad", grad
        maxslope=float(MAXGRAD-grad)/3.
        print "grad ", grad, ", max slope ",maxslope
        for ind2 in range(0,slopesteps+1):
            slope=miny+(maxy-miny)/float(slopesteps)*ind2
            if slope>(maxslope+0.2):
                #print "mask ", ind1, ind2, grad, slope
                masking[ind2, ind1]=1

    #print masking
    masked = np.ma.masked_where(masking == 0, masking)
    return masking

def makemask1(X, Y):
    #calculate max slopes for mask
    
    masking=np.zeros( np.shape(X)) #-MINGRAD
    #print masking
    #print maxmax
    #print "looping..."
    for ix,iy in np.ndindex(X.shape):
        maxslope=float(MAXGRAD-X[ix,iy])/3.
        #print "grad ", grad, ", max slope ",maxslope
        #print ix, iy
        if Y[ix,iy]>maxslope+0.2:
            #print X[ix, iy], Y[ix, iy], maxslope
            masking[ix,iy]=1

    #print masking
    masked = np.ma.masked_where(masking == 0, masking)
    return masked

def makemask2(X, Y,Z):
    #calculate max slopes for mask
    
    masking=np.zeros( np.shape(X)) #-MINGRAD
    #print masking
    #print maxmax
    #print "looping..."
    for ix,iy in np.ndindex(X.shape):
        maxslope=float(MAXGRAD-X[ix,iy])/3.
        #print "grad ", grad, ", max slope ",maxslope
        #print ix, iy
        if Y[ix,iy]>maxslope+0.2:
            #print X[ix, iy], Y[ix, iy], maxslope
            masking[ix,iy]=1

    #print masking
    masked = np.ma.masked_where(masking == 1, Z)
    return masked
#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################

fig, ax = plt.subplots()

if len(sys.argv) <2:
    print "This is the program 'plot_networkpars.py'"
    print "Usage: ./plot_networkpars.py <figure name> <input file names ...>"
    sys.exit(1)
else:
    figname=sys.argv[1]

count=0


#print masked
for filename in sys.argv[2:]:
    
    pars = np.loadtxt(filename, dtype='f', delimiter='\t',ndmin=2)

    #print pars

    sizes = np.shape(pars[:,:])
    print "sizes ", sizes
    if(sizes[0]>1):
        minInColumns = np.amin(pars, axis=0)
        maxInColumns = np.amax(pars, axis=0)
        print minInColumns
        print maxInColumns
        #
    
    
    if not count:
        X, Y, Z = density_estimation(pars[:,0], pars[:,1],minInColumns,maxInColumns)
    
        levels = MaxNLocator(nbins=20).tick_values(Z.min(), Z.max())
        # pick the desired colormap, sensible levels, and define a normalization
        # instance which takes data values and translates those into levels.
        cmap = plt.get_cmap('gist_earth_r')
        cmap2 = plt.get_cmap('binary')
        masked=makemask2(X,Y,Z) #this one masks parts of Z itself, rather than pasting an image on top of the image of Z
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

        #plot the density heat map
        im = ax.pcolormesh(X,Y,masked, cmap=cmap, norm=norm)
        #plot contour lines
        #ax.contour(X, Y, masked, colors='black')             

        #plot the mask  
        #masking=makemask(X,Y)
        #print masking
        #levels2 = MaxNLocator(nbins=20).tick_values(masking.min(), masking.max())
        #norm2=BoundaryNorm(levels2, ncolors=cmap2.N, clip=True)
        #im2=ax.pcolormesh(X,Y,masking,alpha=1, cmap=cmap2,)#norm=norm2
        ax.set_ylim(minInColumns[1],maxInColumns[1])
        ax.set_xlim(minInColumns[0]-50,maxInColumns[0])
 #       ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,extent=[minInColumns[0]-50, maxInColumns[0], minInColumns[1], maxInColumns[1]]) 
        
 #       ax.imshow(masking, cmap='binary',extent=[minInColumns[0]-50, maxInColumns[0], minInColumns[1], maxInColumns[1]], origin="lower")

    ax.scatter(pars[:,0], pars[:,1], s=6+count*10,c=colours[count])   
    count+=1


ax.set_aspect('auto')                                                              
fig.savefig(figname)
plt.show()