#script to read the W AND WA matrix from an elli file and make a dot file with circle layout. 
import sys
from pylab import *
from math import *
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import random as random
from collections import Counter

inprops=["gradient", "divisions"] #, "ccs"] 

outprops=["divisions", "none"] #, "ccs"]

maxreg=1.5

colours=["salmon", "lightskyblue", "palegreen"]

#colours=['r', 'b']
gen1=[] #for storing w matrix
gen2=[] #for stats on outdegree
genlab=[] #for gene labels
############
### Main ###
############

#read data file
#and open write file (for network)
file1=open(sys.argv[1],'r')
file2=open(sys.argv[2],'r')

networkfilename=sys.argv[3]+"_diff.dot"
file3=open(networkfilename,'w')

#write standard stuff to network file
file3.write("digraph G {\n\n")
file3.write("  layout=\"dot\"\n")
file3.write("  node [fontname=\"arial\",fontsize=18,style=filled];\n\n fixedsize=true;\n size=\"10,10\";\n resolution=100;\n bgcolor=\"#FFFFFF\";\n\n")

#first read labels if appropriate
check=0
gen=0

count=0
gencount=1


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
      print "Mismatch between specified and actual nr inputs. proceed with caution..."

    if (outnr1 != len(outprops) or outnr2 != len(outprops)):
      print "Mismatch between specified and actual nr outputs. proceed with caution..."

  elif (count==1):
    for el1,el2 in zip(arr1, arr2):
      file3.write(str(gencount)+"[label=\""+inprops[gencount-1]+" "+str(float(el2)-float(el1))+"\" color=\""+colours[0]+"\"];\n")
      gencount+=1

  else: 
    gennr=int(arr1[1])+1
    if (int(arr1[0])==1): #regulatory node
      file3.write(str(gennr)+"[label=\""+str(gennr)+" "+str(float(arr2[2])-float(arr1[2]))+"\" color=\""+colours[int(arr1[0])]+"\"];\n")
      gencount=1
      for el1,el2 in zip(arr1[3:], arr2[3:]):
        diff=float(el2)-float(el1)
        if diff>0.0000001:
          file3.write(str(gencount)+"->"+str(gennr)+"[label=\""+str(diff)+"\" color=\"royalblue\" penwidth="+str(max(0.5,10.0*diff/maxreg))+"];\n")
        elif diff<-0.0000001:
          file3.write(str(gencount)+"->"+str(gennr)+"[label=\""+str(diff)+"\" color=\"firebrick\" penwidth="+str(max(0.5,10.0*abs(diff)/maxreg))+"];\n")
        gencount+=1
    else:
      file3.write(str(gennr)+"[label=\""+outprops[gennr-1-innr1-regnr1]+" "+arr1[2]+"\" color=\""+colours[int(arr1[0])]+"\"];\n")
      gencount=1+innr1
      for el1,el2 in zip(arr1[3:], arr2[3:]):
        diff=float(el2)-float(el1)
        if diff>0.0000001:
          file3.write(str(gencount)+"->"+str(gennr)+"[label=\""+str(diff)+"\" color=\"royalblue\" penwidth="+str(max(0.5,10.0*diff/maxreg))+"];\n")
        elif diff<-0.0000001:
          file3.write(str(gencount)+"->"+str(gennr)+"[label=\""+str(diff)+"\" color=\"firebrick\" penwidth="+str(max(0.5,10.0*abs(diff)/maxreg))+"];\n")
        gencount+=1
     
    
 
  count+=1

file3.write("}")

