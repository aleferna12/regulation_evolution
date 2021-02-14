
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

colours=["aqua", "red"]
MINGRAD=50
MAXGRAD=1550

xedges = np.arange(-0.5, 29+0.5, 1.) #range is 0-1550
#yedges = np.arange(-1./6., 9.+1./6., 1./3.)
yedges = np.arange(-7.-1/6., 9.+1./6., 1./3.)

# from: https://stackoverflow.com/questions/57277247/the-point-that-minimizes-the-sum-of-euclidean-distances-to-a-set-of-n-points
def geometric_median(pts, numIter = 1200):
    from numpy.linalg import norm as npnorm
    c_pt_old = np.array([np.median(pts[:,0]),np.median(pts[:,1])])
    c_pt_new = np.array([0,0])

    while npnorm(c_pt_old-c_pt_new)>1.e-2:
        num = 0
        denom = 0
        #print c_pt_old, c_pt_new
        for i in range(len(pts)):
            dist = npnorm(c_pt_new-pts[i,:])
            if dist>0.00000000000001:
                num += pts[i,:]/dist
                denom += 1/dist
        c_pt_old = c_pt_new
        c_pt_new = num/denom

    return c_pt_new

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

    if not count:
        H, xedges, yedges = np.histogram2d(pars[:,0], pars[:,1], bins=(xedges, yedges))
        #print xedges
        X, Y = np.meshgrid(xedges, yedges) # <- invert x,y
        #print X
        cmap = plt.get_cmap('gist_earth_r')
        
        ax.pcolormesh(X[:-1], Y[:-1], H.T, cmap = cmap)
        #ax.scatter(pars[:,0],pars[:,1], s=12,c="darkgray")   
        gmedian=geometric_median(pars)
        print gmedian[0], gmedian[1]
        ax.scatter(gmedian[0],gmedian[1], s=16,c=colours[count])   
        
    if count:
        ax.scatter(pars[:,0], pars[:,1], s=16,c=colours[count])   
    count+=1


ax.set_aspect('auto')                                                              
fig.savefig(figname)
plt.show()


