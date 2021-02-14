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
networkfilename=sys.argv[1][:sys.argv[1].rfind(".")]+"_network.dot"
file2=open(networkfilename,'w')

#write standard stuff to network file
file2.write("digraph G {\n\n")
file2.write("  layout=\"dot\"\n")
file2.write("  node [fontname=\"arial\",fontsize=18,style=filled];\n\n fixedsize=true;\n size=\"10,10\";\n resolution=100;\n bgcolor=\"#FFFFFF\";\n\n")

#first read labels if appropriate
check=0
gen=0

count=0
gencount=1


for line in file1:

  arr1=line.split()

  if (not count):
  
    innr=int(arr1[0])
    regnr=int(arr1[1])
    outnr=int(arr1[2])

    if (innr != len(inprops)):
      print "Mismatch between specified and actual nr inputs. proceed with caution..."

    if (outnr != len(outprops)):
      print "Mismatch between specified and actual nr outputs. proceed with caution..."

  elif (count==1):
    for el in arr1:
      file2.write(str(gencount)+"[label=\""+inprops[gencount-1]+" "+el+"\" color=\""+colours[0]+"\"];\n")
      gencount+=1

  else: 
    gennr=int(arr1[1])+1
    if (int(arr1[0])==1): #regulatory node
      file2.write(str(gennr)+"[label=\""+str(gennr)+" "+arr1[2]+"\" color=\""+colours[int(arr1[0])]+"\"];\n")
      gencount=1
      for el in arr1[3:]:
        if float(el)>0.0000001:
          file2.write(str(gencount)+"->"+str(gennr)+"[label=\""+el+"\" color=\"royalblue\" penwidth="+str(max(0.5,10.0*float(el)/maxreg))+"];\n")
        elif float(el)<-0.0000001:
          file2.write(str(gencount)+"->"+str(gennr)+"[label=\""+el+"\" color=\"firebrick\" penwidth="+str(max(0.5,10.0*math.sqrt(float(el)*float(el))/maxreg))+"];\n")
        gencount+=1
    else:
      file2.write(str(gennr)+"[label=\""+outprops[gennr-1-innr-regnr]+" "+arr1[2]+"\" color=\""+colours[int(arr1[0])]+"\"];\n")
      gencount=1+innr
      for el in arr1[3:]:
        if float(el)>0.0000001:
          file2.write(str(gencount)+"->"+str(gennr)+"[label=\""+el+"\" color=\"royalblue\" penwidth="+str(max(0.5,10.0*float(el)/maxreg))+"];\n")
        elif float(el)<-0.0000001:
          file2.write(str(gencount)+"->"+str(gennr)+"[label=\""+el+"\" color=\"firebrick\" penwidth="+str(max(0.5,10.0*math.sqrt(float(el)*float(el))/maxreg))+"];\n")
        gencount+=1
     
    
 
  count+=1

file2.write("}")

