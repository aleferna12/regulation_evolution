#script to read the W AND WA matrix from an elli file and make a dot file with circle layout. 
import sys
from math import *
import numpy as np
from glob import glob
import random as random
from collections import Counter

inprops=["gradient", "divisions"] #, "ccs"] 

outprops=["divisions"] #, "ccs"]

inputrange=[1,1] #[1500, 3]

############
### Main ###
############

#read data file
#and open write file (for network)
file1=open(sys.argv[1],'r')
file2=open(sys.argv[2],'r')

#file3=open(networkfilename,'w')

#first read labels if appropriate
check=0
gen=0

count=0
gencount=1

totaldiff=0.

for line1,line2 in zip(file1,file2):

  arr1=line1.split()
  arr2=line2.split()
  
  if (not count):
  
    innr1=int(arr1[0])
    regnr1=int(arr1[1])
    outnr1=int(arr1[2])

    innr2=int(arr2[0])
    regnr2=int(arr2[1])
    outnr2=int(arr2[2])

    if (innr1 != len(inprops) or innr2 != len(inprops)):
      print ("Mismatch between specified and actual nr inputs. proceed with caution...")

    if (outnr1 != len(outprops) or outnr2 != len(outprops)):
      print ("Mismatch between specified and actual nr outputs. proceed with caution...")

  elif (count==1):
    #print arr1
    #print arr2
    for ind, (el1,el2) in enumerate(zip(arr1, arr2)):
      totaldiff+=((float(el2)-float(el1))**2)*inputrange[ind]
      gencount+=1

  elif count<6: 
    gennr=int(arr1[1])+1
    if (int(arr1[0])==1): #regulatory node
      totaldiff+=(float(arr2[2])-float(arr1[2]))**2
      gencount=1
      for el1,el2 in zip(arr1[3:], arr2[3:]):
        totaldiff+=(float(el2)-float(el1))**2
        gencount+=1
    else:
      gencount=1+innr1
      for el1,el2 in zip(arr1[3:], arr2[3:]):
        totaldiff+=(float(el2)-float(el1))**2
        gencount+=1
  count+=1

print ("{}".format(totaldiff**0.5))
#file3.write(str(totaldiff))
#file3.close()

