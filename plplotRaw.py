"""
  Function:  plplotRaw
  --------------------
  Plot the power-law function, take raw (not cumulative) occurences.
  
  momo.sander@ebi.ac.uk
"""     

import matplotlib.pyplot as plt
from math import *

def plplotRaw(x, filename):

  filename = filename.split('\t')[0]
  y = []
  xu = unique(x)
  xBin = []
  for valu in xu:
    valb = 0
    for val in x:
      if  valu == val:
        valb +=1
    xBin.append(valb)

  c1 = xu
  c2 = xBin
  F = plt.figure(1)
  h =plt.loglog(c1, c2, 'bo',marker = 'x', markersize=8,markerfacecolor=[1,1,1],markeredgecolor=[0,0,1])

  xr1 = pow(10,floor(log(min(x),10)))
  xr2 = pow(10,ceil(log(min(x),10)))
  yr2 = max(xBin)
 

  plt.axhspan(ymin=0.5,ymax=yr2,xmin=xr1,xmax=xr2)
  plt.ylabel('Pr(X >= x)',fontsize=12);
  plt.xlabel('x',fontsize=12)
  plt.subplots_adjust(left=0.3, bottom=0.3)
  F.set_size_inches(3,2)
  F.savefig('visual/plplot%s.pdf' % filename)
  plt.clf()               
  return h

# helper function unique
def unique(seq): 
    # not order preserving 
    set = {} 
    map(set.__setitem__, seq, []) 
    return set.keys()



