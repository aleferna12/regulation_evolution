#!/usr/bin/python2.7

'''
take a key lock file and returns jvalues
'''

import sys,math,os,subprocess
#from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np

def JWithMedium(key):
  Jval=0
  for i in range(keypos_formedium):
    Jval += int(key[i])*pow(2,keypos_formedium-i-1); #so that zeroth bit is most significant
  Jval += 10; #so that interaction with medium can't be 0
  return Jval

def JWithOtherTau(( key1,lock1 ),( key2,lock2 )):
  score=0;
  
  for i in range(len(key1)):
    score += 1 if key1[i] != lock2[i] else 0;
    score += 1 if key2[i] != lock1[i] else 0;
  
  Jval = 3 + (int)(0.5+ 40.*math.exp( -0.01 * pow(float(score),2.) ));
  return Jval
  
#filename="data_cellcount.txt"

keypos_formedium=4

if len(sys.argv)==1:
  print "need a file name with -filename option"
  sys.exit(1)
elif len(sys.argv)>1:
  pos=1
  while pos < len(sys.argv) :
    #if sys.argv[pos]=='-maxtime':
    #  pos+=1
    #  maxTime=float(sys.argv[pos])
    if sys.argv[pos]=='-filename':
      print "Changing filename to",
      pos+=1
      filename=sys.argv[pos]
      print filename
    else:
      print "This is the program 'plot_cellcount.py'"
      print "It plots files formatted as time countpreys countpredators"
      print "Here are possible options:"
      print "-maxtime [INT]"
      print "-filename [filename]"
      sys.exit(1)
    pos +=1

d1={}
with open(filename,"r") as fin:
  lines = fin.readlines()
  taus=int(lines[0])
  howmany=len(lines)
  i=1
  tau=1
  while i < howmany:
    key = "".join(lines[i].split()) 
    i+=1
    lock= "".join(lines[i].split())
    i+=1
    d1[tau]=(key,lock)
    tau+=1

print taus
print d1
print 

for key in sorted(d1):
  print "Tau =",key
  print "Half J with medium:",JWithMedium(d1[key][0])
  for key2 in sorted(d1):
    print "J with tau:",key2,"=",JWithOtherTau(d1[key],d1[key2])
  print  


