import sys
import matplotlib.pyplot as plt
import numpy as np
import os

# this program reads input from a script which has assessed how networks react to a particular combination of gradient and division status
# the script has produced for each network a matrix with 0 (migrate) and 1 (divide), which this program will plot and find the consensus for.


if len(sys.argv) <3:
    print ("This is the program 'plot_netanalysis_jan.py'")
    print ("Usage: ./plot_netanalysis_jan.py <output_file> <plot individuals?> <input filenames>")
    sys.exit(1)
else:
    outputfile=sys.argv[1]
    indiplot=int(sys.argv[2])

arraystorage=[]
filestorage=[]
init=0
count=0
sizes=None
consensus=None

for filename in sys.argv[3:]:
    #print ("{}".format(filename))
    divmig = np.loadtxt(filename, dtype='i', delimiter='\t')
       
    #print sizes 
    if not init:
        sizes = np.shape(divmig[1:,1:])
        consensus=np.zeros((sizes[0]*sizes[1],),dtype=int)
        init=1

    outfile=os.path.splitext(filename)[0]
    #for if you still need to plot the individuals:
    if (indiplot):
        fig=plt.figure() #!
        fig.set_size_inches(1, 1.*sizes[0]/sizes[1], forward = False) #!
        ax = plt.Axes(fig, [0., 0., 1., 1.]) #!
        ax.set_axis_off() #!
        fig.add_axes(ax) #!
        ax.imshow(divmig[1:,1:], cmap='RdYlBu', origin='lower')
    
        divshare=divmig[1:,1:].sum()
        migshare=(sizes[0])*(sizes[1])-divshare
        
        migs="%04d" % (migshare,)
        #print divs
        plt.savefig("div_"+str(migs)+"_"+outfile+".pdf", dpi=sizes[1]) #bbox_inches='tight'
        plt.close()

    binarystring=divmig[1:,1:].flatten()

    consensus=np.add(binarystring, consensus)
    #print ("{}".format(consensus))
    arraystorage.append(binarystring)
    filestorage.append(outfile)

    count+=1




#find the consensus sequence
bool_consensus= consensus > count/2
print ("{}".format(bool_consensus))
consensus_sequence=bool_consensus.astype(int)
print ("consensus is {}".format(consensus_sequence))
wfilename="consensussequence_"+outputfile+".dat"
writefile=open(wfilename,"w")
for el in consensus_sequence:
    writefile.write(str(el)+" ")
writefile.close()

#display consensus image
imcons=np.reshape(consensus_sequence,sizes)
fig=plt.figure() #!
fig.set_size_inches(1, 1.*sizes[0]/sizes[1], forward = False) #!
ax = plt.Axes(fig, [0., 0., 1., 1.]) #!
ax.set_axis_off() #!
fig.add_axes(ax) #!
ax.imshow(imcons, cmap='RdYlBu', origin='lower')

#outfile=os.path.splitext(outputfile)[0]
plt.savefig("consensus"+"_"+outputfile+".pdf", dpi=sizes[1]) #bbox_inches='tight'
plt.close()

#find for each individual the distance to the consensus sequence
#writefile=open(outputfile, "w")
#fig=plt.figure() #
#hamms=[]
minhamm=999999999
for fi,seq in zip(filestorage, arraystorage):
    hamm=np.count_nonzero(seq!=consensus_sequence)
    if hamm<minhamm:
        minhamm=hamm
        minfile=fi

print ("file with individual closest to consensus: {}".format(minfile))
#    hamms.append[hamm]
    #writefile.write(fi+"\t"+str(hamm)+"\n")
#maxbina=max(hamms)
#hista, bin_edgesa = np.histogram(hamms, bins = range(maxbina))
#plt.plot(bin_edgesa[:-1],hista)

#writefile.close()




