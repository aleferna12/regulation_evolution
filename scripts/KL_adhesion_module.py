
'''
contains standard adhesion module for cell evolution type experiments.
In the future it should also be able to read par files
'''

import sys,math,os,subprocess
import numpy as np

lenkl=24

MIN_Jval_withmedium = 8
lookuptable_Jmedium=[4,3,2,1,1,1]
keypos_formedium=len(lookuptable_Jmedium)

def Distance(current,target):
  return sum([ abs(x-y) for x,y in zip(target,current) ])

# in this new version we just use the lookup table    
def ReverseKeyForMedium(hjm,lenkl,lookuptable_Jmedium):
  if hjm<MIN_Jval_withmedium or hjm>MIN_Jval_withmedium+sum(lookuptable_Jmedium):
    print ("Error, value for J(cell,medium) ",hjm," is out of range [",MIN_Jval_withmedium,MIN_Jval_withmedium+sum(lookuptable_Jmedium),"]")
    sys.exit(1)
  
  key=lenkl*[0]
  hjm-=10
  i=keypos_formedium-1
  while i>=0:
    #print i,key
    key[i]=hjm%2
    hjm=hjm/2
    i-=1
  return key

def JWithMedium(key,lookuptable_Jmedium):
  #Jval=0
  #for i,j in xrange(len(lookuptable_Jmedium)):
    #Jval += int(key[i])*pow(2,keypos_formedium-i-1); #so that zeroth bit is most significant
  Jval = sum( [ x*y  for x,y in  zip( key[:len(lookuptable_Jmedium)] , lookuptable_Jmedium  ) ]  )
  Jval += MIN_Jval_withmedium; #so that interaction with medium can't be 0
  return Jval

####   GENERALISED FORMULA   ####
#def JWithOtherTau( ( key1,lock1 ),( key2,lock2 ) ):
def JWithOtherTau( tup1,tup2 ):
  score=0;
  key1=tup1[0]
  lock1=tup1[1]
  key2=tup2[0]
  lock2=tup2[1]
    
  for i in range(len(key1)):
    score += 1 if key1[i] != lock2[i] else 0;
    score += 1 if key2[i] != lock1[i] else 0;
  
  #Jval = 3 + (int)(0.5+ 40.*math.exp( -pow(float(score),2.) / pow( float(lenkl) , 2.) ));
  #the above formula has been changed to 
  Jval =(int)(0.5+  52.- 48./float(2*lenkl) *float(score))
  #print("Jval: ", Jval)
  return Jval
  
def RandomKL(lenkl):
  rk1=[ 1 if np.random.rand()<0.5 else 0 for _ in range(lenkl) ] # randomises everything
  rk1=[ 1 if np.random.rand()<0.5 else 0 for _ in range(lenkl) ]
  
  return rk1,rl1

def RandomiseKL(k1,l1,k2,l2,mutrate):

  rk1=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(k1) ] # randomises everything
  rk2=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(k2) ] # randomises everything
  #rk1=[ 1-x if (np.random.rand()<mutrate and i>=keypos_formedium) else x for i,x in enumerate(k1) ]
  #rk2=[ 1-x if (np.random.rand()<mutrate and i>=keypos_formedium) else x for i,x in enumerate(k2) ]
  
  rl1=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(l1) ]
  rl2=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(l2) ]
  #print (rk1, rl1)
  return rk1,rl1,rk2,rl2
