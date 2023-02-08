#!/usr/bin/python2.7

'''
Generates a bunch of KL and calculate adhesion
'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np

# lenkl=25
# howmany = 10
# alphabet

# [5,4,3,1,1,1] -> means [0->15]
# [4,3,2,1,1,1] -> means [0->12]
def JWithMedium(key,lookuptable_Jmedium):
  return 20
  #Jval=0
  #for i,j in xrange(len(lookuptable_Jmedium)):
    #Jval += int(key[i])*pow(2,keypos_formedium-i-1); #so that zeroth bit is most significant
  
  # Jval = sum( key[:6]  )
  # Jval += 2
  # return Jval
  # 
  newlist = [ 1 if x>1 else 0 for x in key[:len(lookuptable_Jmedium)]]
  Jval = sum( [ x*y  for x,y in  zip( newlist , lookuptable_Jmedium  ) ]  )
  Jval += 8
  return Jval


####   GENERALISED FORMULA   ####
def JWithOtherTau(( key1,lock1 ),( key2,lock2 )):
  score=0;
    
  # for i in range(len(key1)):
  #   score += 1./12.*(i+1) if key1[i] != lock2[i] else 0;       #1./12.*i
  #   score += 1./12.*(i+1) if key2[i] != lock1[i] else 0;
  # 
  
  bestscore=0
  
   
  for i in range(len(key1)):
    score += 1. if key1[i] == lock2[i] else 0;       #1./12.*i
    score += 1. if key2[i] == lock1[i] else 0;
  
    #score += 1. if key1[i] != lock2[i] else 0;       #1./12.*i
    #score += 1. if key2[i] != lock1[i] else 0;
  
  #Jval =int( 52.- 48./float(2*lenkl) *float(score))
  #Jval =int( 52.- float(score))
  Jval=72-score
  return Jval
   
  
  # Jval = 5 + (int)(0.5+ 40.*math.exp( -pow(float(score),2.) / pow( float(lenkl) , 2.) ));
  # print score
  # return score
  
  # x=float(score)
  # slope=0.3
  # Jval =44.+12.-  (12.*math.exp(lenkl*slope) + 44.*math.exp(slope*x))/(math.exp(lenkl*slope) + math.exp(slope*x))
  # 
  
  # Jval =4.+24.+24.*math.cos(math.pi*float(score/4.)/48.)
  # Jval = 4.+ (2.*lenkl)/(1.+ (float(score)/float(lenkl))**2.)
  #Jval = 4+ (28-score)
  
  return Jval
  
#filename="data_cellcount.txt"

def RandomiseKL(k1,l1,mutrate):
  alphabet = [0,1]
  # alphabet = [0,1,2,3]
  #rk1=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(k1) ] # randomises everything
  #rl1=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(l1) ] # randomises everything
  
  #rk1=[ random.choice(alphabet[ : alphabet.index(x) ]+  alphabet[ 1+alphabet.index(x): ] ) if np.random.rand()<mutrate else x for i,x in enumerate(k1) ] # randomises everything
  #rl1=[ random.choice(alphabet[ : alphabet.index(x) ]+  alphabet[ 1+alphabet.index(x): ] ) if np.random.rand()<mutrate else x for i,x in enumerate(l1) ] # randomises everything
  # print k1
  
  
  rk1=[ random.choice(alphabet) for i,x in enumerate(k1) ] # randomises everything
  rl1=[ random.choice(alphabet) for i,x in enumerate(l1) ] # randomises everything
  # print rk1
  #rk1=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(k1) ] # randomises everything
  
  #rl1=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(l1) ]
  #rl2=[ 1-x if np.random.rand()<mutrate else x for i,x in enumerate(l2) ]
  return rk1,rl1

#######################################
####### BEGIN STUFF ###################
#######################################

#keypos_formedium=4
# lenkl=96
#lenkl=24
lenkl=24
N=100000

k1=lenkl*[0]
k2=lenkl*[0]
l1=lenkl*[0]
l2=lenkl*[0]
mutrate = 0.5

#lookuptable_Jmedium=[5,3,2,1,1,1,1,1]
#lookuptable_Jmedium=[5,4,3,1,1,1]
lookuptable_Jmedium=[4,3,2,1,1,1]
keypos_formedium=len(lookuptable_Jmedium)

lG11=[]
lG12=[]
lG21=[]
lG22=[]
lJ10=[]
lJ11=[]
lJ12=[]

for i in xrange(N):
    if(i%10000==0): print i
    k1,l1 = RandomiseKL(k1,l1,mutrate)
    k2,l2 = RandomiseKL(k2,l2,mutrate)
    
    J10 = JWithMedium(k1,lookuptable_Jmedium)
    J20 = JWithMedium(k2,lookuptable_Jmedium)
    J12 = JWithOtherTau(( k1,l1 ),( k2,l2 ))
    J11 = JWithOtherTau(( k1,l1 ),( k1,l1 ))
    J22 = JWithOtherTau(( k2,l2 ),( k2,l2 ))
    
    lJ10.append(J10)
    lJ11.append(J11)
    lJ12.append(J12)
    
    lG11.append (J10 - J11/2.)
    #lG22.append (J20 - J22/2.)
    
    lG12.append (J10 - J12/2.) 
    #lG21.append( J20 - J12/2.)
    
f, (ax1, ax2, ax3, ax4,ax5) = plt.subplots(5, 1, sharex=True)
bins = 100
#bins = np.arange( min(lG11)- 1., max(lG11) + 1.5) - 0.5
#ax1.hist(lG11,bins, label="$\gamma 11$") 
ax1.scatter(lG11,[ 100.*random.random() for x in lG11 ] ,s=0.1, alpha=0.1, label="$\gamma 11$")
ax1.legend()

#bins = np.arange( min(lG12)- 1., max(lG12) + 1.5) - 0.5
# ax2.hist(lG12, bins, label="$\gamma 12$")    
ax2.scatter(lG12,[ 100.*random.random() for x in lG12 ] ,s=0.1, alpha=0.1, label="$\gamma 12$")

ax2.legend()

#bins = np.arange(0, max(lJ10) + 1.5) - 0.5
# ax3.hist(lJ10,bins, label="J10")
ax3.scatter(lJ10,[ 100.*random.random() for x in lJ10 ] ,s=0.1, alpha=0.1, label="J10")

ax3.legend()

#bins = np.arange(0, max(lJ11) + 1.5) - 0.5
# ax4.hist(lJ11, bins, label="J11")    
ax4.scatter(lJ11,[ 100.*random.random() for x in lJ11 ] ,s=0.1, alpha=0.1, label="J11")

ax4.legend()

#bins = np.arange(0, max(lJ12) + 1.5) - 0.5
# ax5.hist(lJ12, bins, label="J12")
ax5.scatter(lJ12,[ 100.*random.random() for x in lJ12 ] ,s=0.1, alpha=0.1, label="J12")

ax5.legend()

plt.show()

