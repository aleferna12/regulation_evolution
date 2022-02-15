##!/usr/bin/python2.7

'''
take jvalues and return key lock file only for ONE cell type

Here I try to do it with length 20
'''

import sys,math,os,subprocess
#from PIL import Image
import KL_adhesion_module as klam

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np

def Distance(current,target):
  return sum([ abs(x-y) for x,y in zip(target,current) ])

def PrintSuccess(target,current,k1,l1):
  print ("")
  print ("Success")
  PrintResult(target,current,k1,l1)
  
def Print5000(target,current,k1,l1):
  print ("")
  print (howmanysteps,"steps passed, this is what I got")
  PrintResult(target,current,k1,l1)
  
def PrintResult(target,current,k1,l1):
  print ("targ:", target)
  print ("curr:",current)
  print ("distance", Distance(current,target))
  print ("with following K1,L1 (for KL...dat file)")
  print ("")
  print ("3")
  for x in k1:
    print (x,end = '')
  print ("") 
  for x in l1:
    print (x,end = '')
  print ("")


def RandomiseKL(k1,l1,mutrate):

  rk1=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(k1) ] # randomises everything
  #rk2=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(k2) ] # randomises everything
  
  rl1=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(l1) ]
  #rl2=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(l2) ]
  return rk1,rl1

#######################################
####### BEGIN STUFF ###################
#######################################

#keypos_formedium=4
# klam.lenkl
# klam.lookuptable_Jmedium

keypos_formedium=len(klam.lookuptable_Jmedium)

argi=0
hj1m = int( sys.argv[1]) #10 #up to 15
j11 = int( sys.argv[2]) #6 # to 41 #int(sys.argv[1])  # variable of interest

#j12 = int( sys.argv[3]) #8 # to 14  #high adh with prey - always
#hj2m = hj1m
#j21 = j12
#j22 = int( sys.argv[4]) #38 # up to 41

#while j11<43:
    #sys.stderr.write("j11="+str(j11)+"\n")
    
    
    #hj1m = 10 #up to 15
    #j12 = 8 # to 14  #high adh with prey - always
    #hj2m = hj1m
    #j21 = j12
    #j22 = 38 # up to 41

mutrate=0.0500
howmanysteps=50000

k1=klam.lenkl*[0]
#k2=lenkl*[0]

k1=klam.ReverseKeyForMedium(hj1m,klam.lenkl,klam.lookuptable_Jmedium)
#k2=ReverseKeyForMedium(hj2m,lenkl,lookuptable_Jmedium)

k1=[ int(2.*np.random.rand()) if i>=keypos_formedium else x for i,x in enumerate(k1) ]
#k2=[ int(2.*np.random.rand()) if i>=keypos_formedium else x for i,x in enumerate(k2) ]

l1=[ int(2.*np.random.rand()) for _ in range(klam.lenkl) ]
#l2=[ int(2.*np.random.rand()) for _ in range(lenkl) ]

#print k1,k2,l1,l2
#exit(1)

#print k1
#print k2
#print l1
#print l2

#Now the cool part:
# calculate product with vector made by (hj1m,j11 ,j12 ,hj2m ,j21 ,j22 )
# if it is very close, ok, else randomly change something, etc...

# no
#the max possible score is sum(hj1m^2,j11^2 ,j12^2 , etc...)
#maxval= sum( [hj1m**2,j11**2 ,j12**2 ,hj2m **2,j21**2 ,j22**2] )

#init_targetJ = [hj1m,j11 ,j12 ,hj2m ,j21 ,j22 ]
init_targetJ = [hj1m,j11 ]
print ("Target: ", [hj1m,j11 ], "pos 0 -> +5, pos 2,4 -> +6, pos 5 -> +3")
currentJ = [ klam.JWithMedium(k1,klam.lookuptable_Jmedium), 
             klam.JWithOtherTau(( k1,l1 ),( k1,l1 ))
           ]
  
distance = klam.Distance(init_targetJ,currentJ)

counter=0
counter2=0

targetJ=init_targetJ[:]

lresults=[]

          
counter=0

#print "Target:",targetJ

distance=klam.Distance(targetJ,currentJ) 

while True:
          if counter%10000==0: print (counter) #,end = '')
          if distance==0:
            lresults.append( [ distance, currentJ, k1,l1] )
            break
            
          if counter>=howmanysteps:
            #Print5000(targetJ,currentJ,k1,l1,k2,l2)
            lresults.append( [ distance, currentJ, k1,l1] )
            print (lresults[-1])
            #print "Got :", distance, currentJ
            break
          
          counter+=1
          #Randomise 
          rk1,rl1 = RandomiseKL(k1,l1,mutrate)
          while (rk1,rl1) == (k1,l1):
            rk1,rl1 = RandomiseKL(k1,l1,mutrate)
          rcurrentJ = [ klam.JWithMedium(rk1,klam.lookuptable_Jmedium), 
                      klam.JWithOtherTau(( rk1,rl1 ),( rk1,rl1 ))
                      ]
          
          rdistance = klam.Distance(targetJ,rcurrentJ) 
          
          if rdistance<=distance:
            k1=rk1[:]
            l1=rl1[:]
            print (k1, l1)
            distance=rdistance
            currentJ=rcurrentJ[:]
            
            
print ("I'm here!")
best=[[100000,'bla']]
for result in lresults:
  if result[0]<best[0][0]:
    best=[result[:]]
  elif result[0]==best[0][0]:
    best.append( result[:] )
    
print ("")
for bla in best:
  print (bla[0])
  print (bla[1])
  print ("here2")
  for x in bla[2:]:
    for y in x: print (y,end = '')
    print ("")
  print ("")

j11+=1
