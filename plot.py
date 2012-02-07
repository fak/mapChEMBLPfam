
"""
Function:  plotHist 

from a List containing the raw observations (eg. x= [1,0,1,2,3,2,1,1,2,1,1,2,1]
plot a basic histogram directly from python.
  --------------------

"""

def plotHist(x,y,app):
  import matplotlib.pyplot as plt
  import numpy as np
  
  plt.grid(True)  
  plt.hist([x,y], max(x), normed = True, label = ['genome', 'ChEMBL'])
  plt.xlabel('')
  plt.ylabel('')
  plt.title('')
  plt.xlim(0,12)
  plt.ylim(0,0.6)
  plt.legend()
  plt.savefig('visual/histogram_%s.pdf'%app)
  plt.clf()
  hist, edges = np.histogram(x, bins = [0,1,2,3,4,5,6,7,8,9,10,11,12],normed = True)
  print 'histogram\t:', hist
  print 'edges:\t',edges
  hist, edges = np.histogram(y, bins = [0,1,2,3,4,5,6,7,8,9,10,11,12], normed = True)
  print 'histogram\t:', hist
  print 'edges:\t',edges  
  
def plotEmpCDF(x,app):
  import matplotlib.pyplot as plt
  
  plt.grid(True)  
  plt.hist(x, bins =50, cumulative = True, histtype = 'step' )
  plt.xlabel('residues within predicted domain')
  plt.ylabel('')
  plt.title('')
  plt.xlim(0,1)
  plt.legend()
  plt.savefig('visual/empiricalCDF_%s.pdf'%app)
  plt.clf()
  

def ratioPlots(x,app):
  import matplotlib.pyplot as plt
 
  plt.grid(True)  
  plt.boxplot(x)
  plt.xlabel('')
  plt.ylabel('')
  plt.ylim(0,1.2)
  plt.title('')
  plt.legend()
  plt.savefig('visual/%s.pdf'%app)
  plt.clf()


"""
  Function:  rankPlot
  --------------------
  For each domain, plot rank of occurences in the human genome vs rank in number
  of ligands.
  
  momo.sander@ebi.ac.uk
"""     
def rankPlot(genRankL, ligRankL, rectBords):
  import matplotlib.pyplot as plt
  import matplotlib.patches as ptch
  import numpy as np

    
  xr1 = np.ceil(max(genRankL)/10) 
  yr1 = np.ceil(max(ligRankL)/10)
  yr2 = rectBords  

  F = plt.figure(1)
  plt.plot(genRankL, ligRankL, 'o',markersize=6,markerfacecolor=[0,0,1])
  plt.axhspan(ymin=1,ymax= yr1,xmin=1,xmax=xr1)
  plt.ylabel('Rank in number of ligands',fontsize=12);
  plt.xlabel('Rank in occurences in the human genome',fontsize=12)
  plt.subplots_adjust(left=0.3, bottom=0.3)
  for yr in yr2:
    plt.axhspan(0,yr, facecolor=[1,0,0], alpha = .2)
     
  F.set_size_inches(5,5)
  F.savefig('visual/rankScatter.pdf')
  plt.clf()

  return
  



    

