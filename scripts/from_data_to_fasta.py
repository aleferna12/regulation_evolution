#!/usr/bin/python2.7

import sys,math,os,string,random
from subprocess import Popen, PIPE

"""
takes cell evol formatted files and return fasta files of a certain time step

argument should be file and time step 
"""

#### ARGUMENTS ####

if len(sys.argv)!=3:
  print "Specify some path/to/filename.txt and the time step"
  sys.exit(1)
else:
  for arg in sys.argv[1:]:
    try:
      whichtime=int(arg)
    except:
      filepath=arg
    
print "Time we fetch:", whichtime
print "From file:", filepath

#open filename once to see what is the first time step... for now it is typically zero but in the future?
ldata=[]
success=False
with open(filepath,"r") as fin:
  for line in fin:
    line=line.split()
    time=int(line[0])
    if time<whichtime: continue
    if time>whichtime: 
      if success==False:
        print "time", whichtime,"not found in file", filepath
      break
    if time==whichtime:
      if success!=True:
        success=True;
      ldata.append(line)
    
random.shuffle(ldata)

#we randomise the order of the list
outfile=".".join(filepath.split('.')[:-1]) +"_fasta_time"+str(whichtime)+".fa"
with open(outfile,"w") as fout:
  for i,line in enumerate(ldata):
    fout.write(">"+str(i)+"_time"+str(whichtime)+"_tau"+str(line[1])+"\n")
    bla="".join([ 'A' if x=='0' else 'G' for x in line[3]+line[4] ])
    fout.write(bla+"\n")
    
print "File generated:", outfile
  
