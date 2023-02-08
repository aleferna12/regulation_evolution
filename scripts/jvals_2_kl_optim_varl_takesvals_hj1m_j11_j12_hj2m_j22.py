#!/usr/bin/python2.7

'''
take a key lock file and returns jvalues

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

def PrintSuccess(target,current,k1,l1,k2,l2):
  print
  print "Success"
  PrintResult(target,current,k1,l1,k2,l2)
  
def Print5000(target,current,k1,l1,k2,l2):
  print
  print howmanysteps,"steps passed, this is what I got"
  PrintResult(target,current,k1,l1,k2,l2)
  
def PrintResult(target,current,k1,l1,k2,l2):
  print "targ:", target
  print "curr:",current
  print "distance", klam.Distance(current,target)
  print "with following K1,L1,K2,L2 (for KL...dat file)"
  print 
  print '3'
  for x in k1:
    print x,
  print  
  for x in l1:
    print x,
  print 
  for x in k2:
    print x,
  print  
  for x in l2:
    print x,
  print 


#def ReverseKeyForMedium(hjm,lenkl,lookuptable_Jmedium):
  ##reverse half jvals with medium
  ## inverse formula is: Binary(hj1m - 10)|_4chars
  ## so if we start from half j val 15 -> 5 -> 0101
  ## of course no val higher than 10+Ten(1111)= 25
  ## and min 10+Ten(0000) = 10
  #if hjm<10 or hjm>25:
    #print "hjm out of range [10,25]"
    #sys.exit(1)
  
  #key=lenkl*[0]
  #hjm-=10
  #i=keypos_formedium-1
  #while i>=0:
    ##print i,key
    #key[i]=hjm%2
    #hjm=hjm/2
    #i-=1
  #return key

#######################################
####### BEGIN STUFF ###################
#######################################

if len(sys.argv) != 6:
    print 
    print "*******************************************************************************"
    print " This is  ***   jvals_2_kl_optim_varl_takesvals_hj1m_j11_j12_hj2m_j22.py   *** "
    print "*******************************************************************************"
    print 
    print "Arguments: J(cell1,medium) J(cell1,cell1) J(cell1,cell2) J(cell2,medium) J(cell2,cell2)"
    print "All arguments are necessary, bye!"
    sys.exit(1)

argi=0
hj1m = int( sys.argv[1]) #10 #up to 15
j11 = int( sys.argv[2]) #6 # to 41 #int(sys.argv[1])  # variable of interest
j12 = int( sys.argv[3]) #8 # to 14  #high adh with prey - always
hj2m = int( sys.argv[4])
j21 = j12
j22 = int( sys.argv[5]) #38 # up to 41

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
k2=klam.lenkl*[0]

k1=klam.ReverseKeyForMedium(hj1m,klam.lenkl,klam.lookuptable_Jmedium)
k2=klam.ReverseKeyForMedium(hj2m,klam.lenkl,klam.lookuptable_Jmedium)

k1=[ int(2.*np.random.rand()) if i>=klam.keypos_formedium else x for i,x in enumerate(k1) ]
k2=[ int(2.*np.random.rand()) if i>=klam.keypos_formedium else x for i,x in enumerate(k2) ]

l1=[ int(2.*np.random.rand()) for _ in range(klam.lenkl) ]
l2=[ int(2.*np.random.rand()) for _ in range(klam.lenkl) ]

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

init_targetJ = [hj1m,j11 ,j12 ,hj2m ,j21 ,j22 ]
print "Target: ", [hj1m,j11 ,j12 ,hj2m ,j21 ,j22 ], "pos 0 -> +5, pos 2,4 -> +6, pos 5 -> +3"
currentJ = [ klam.JWithMedium(k1,klam.lookuptable_Jmedium), 
          klam.JWithOtherTau(( k1,l1 ),( k1,l1 )), 
          klam.JWithOtherTau(( k1,l1 ),( k2,l2 )), 
          klam.JWithMedium(k2,klam.lookuptable_Jmedium),
          klam.JWithOtherTau(( k1,l1 ),( k2,l2 )),
          klam.JWithOtherTau(( k2,l2 ),( k2,l2 ))
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
          if counter%10000==0: print counter,
          if distance==0:
            lresults.append( [ distance, currentJ, k1,l1,k2,l2] )
            break
            
          if counter>=howmanysteps:
            #Print5000(targetJ,currentJ,k1,l1,k2,l2)
            lresults.append( [ distance, currentJ, k1,l1,k2,l2] )
            #print "Got :", distance, currentJ
            break
          
          counter+=1
          #Randomise 
          rk1,rl1,rk2,rl2 = klam.RandomiseKL(k1,l1,k2,l2,mutrate)
          while (rk1,rl1,rk2,rl2) == (k1,l1,k2,l2):
            rk1,rl1,rk2,rl2 = klam.RandomiseKL(k1,l1,k2,l2,mutrate)
          rcurrentJ = [ klam.JWithMedium(rk1,klam.lookuptable_Jmedium), 
                      klam.JWithOtherTau(( rk1,rl1 ),( rk1,rl1 )), 
                      klam.JWithOtherTau(( rk1,rl1 ),( rk2,rl2 )), 
                      klam.JWithMedium(rk2,klam.lookuptable_Jmedium),
                      klam.JWithOtherTau(( rk1,rl1 ),( rk2,rl2 )),
                      klam.JWithOtherTau(( rk2,rl2 ),( rk2,rl2 ))
                    ]
          
          rdistance = klam.Distance(targetJ,rcurrentJ) 
          
          if rdistance<=distance:
            k1=rk1[:]
            k2=rk2[:]
            l1=rl1[:]
            l2=rl2[:]
            distance=rdistance
            currentJ=rcurrentJ[:]
            
            
print
best=[[100000,'bla']]
for result in lresults:
  if result[0]<best[0][0]:
    best=[result[:]]
  elif result[0]==best[0][0]:
    best.append( result[:] )
    

for bla in best:
  print bla[0]
  print bla[1]
  for x in bla[2:]:
    for y in x: print y,
    print
  print

j11+=1
