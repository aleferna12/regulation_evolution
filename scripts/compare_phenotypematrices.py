import sys
import matplotlib.pyplot as plt
import numpy as np
import os

#this program compares the phenotype of the individual with which the simulation was started, to (a) matrix/ces from the final population
#this could be the consensus or a representative individual

if len(sys.argv) <1:
    print ("This is the program 'compare_phenotypematrices.py'")
    print ("Usage: ./compare_phenotypematrices.py orifile <input filenames>")
    sys.exit(1)
else:
    orifile=sys.argv[1]

arraystorage=[]
filestorage=[]
init=0
count=0
sizes=None
consensus=None
#print ("starting")
startseq=np.loadtxt(orifile, dtype='i', delimiter=' ')
#print ("{}".format(startseq))

for filename in sys.argv[2:]:
    #print ("{}".format(filename))
    seq = np.loadtxt(filename, dtype='i', delimiter=' ')
    #print   ("{}".format(seq)) 
    hamm=np.count_nonzero(seq!=startseq)
    
    print ("{}".format(hamm))







