"""
  Function:  pdbe
  --------------------
  carries out validation of the algorithm against PDBe.
  
  momo.sander@ebi.ac.uk
"""                                       
def pdbe(pdbDict, release): 
                           

  import numpy as np

  statDict={}
  for target in pdbDict.keys():
    Ts = 0
    Fs = 0
    for cmpdId in pdbDict[target].keys():
      try:
        pdbDict[target][cmpdId]['prediction']
      except KeyError:
        print 'no predictions made for:', target, cmpdId
        continue        
      for pred in pdbDict[target][cmpdId]['prediction']:
        if pred:
          Ts +=1
        else:
          Fs +=1
    tot = Ts+Fs
    if tot > 0:
      statDict[target] = np.true_divide(Ts,tot)  

    
  counts = statDict.values()

  ## Count number of pdb structures in the analysis. 
  pdbCount = {}
  uniprotCount = {}
  ligandCount = {}
  for target in pdbDict.keys():
    uniprotCount[target] = 0
    for ligand in pdbDict[target].keys():
      ligandCount[ligand] = 0 
      for pdb in pdbDict[target][ligand]['pdb']:
        pdbCount[pdb] = 0
  print 'Counts of Uniprot-Ids, molregnos, PDBs',len(uniprotCount.keys()), \
      len(ligandCount.keys()), len(pdbCount.keys()) 
      
  ## Print the histogram.
  hist, edges = np.histogram(counts, bins = 10)
  print 'histogram\t:', hist
  print 'edges:\t',edges
  return counts


"""
  Function:  prepPlot
  --------------------
  Carries out validation of the algorithm against Uniprot.
  
  momo.sander@ebi.ac.uk
"""   
def uniprot(bsDict, release): 
  import numpy as np
  statDict={}
  for target in bsDict.keys():
    Ts = 0
    Fs = 0
    try:
      bsDict[target]['prediction']
    except KeyError:
      print 'no observations made for:', target
      continue        
    for pred in bsDict[target]['prediction']:
      if pred:
        Ts +=1
      else:
        Fs +=1
    tot = Ts+Fs
    if  tot > 0:
      statDict[target] = np.true_divide(Ts,tot)  

  x = statDict.values()
  return x



  """
  Function:  prepPlot
  --------------------
  Prepare to plot a stacked bar to show the True/False mappings in each predList.
  
  momo.sander@ebi.ac.uk
"""     
def prepPlot(predLs, mapTypes):
  import numpy as np
  within = 0.5
  specArr = np.zeros([len(mapTypes), 2])
  for i, mapType in enumerate(mapTypes):
    for predList in predLs[mapType]:
      trues = len(filter(None, predList))
      if np.true_divide(trues, len(predList)) >= within:
        specArr[i][0] += 1
      else:
        specArr[i][1] += 1
        print predList[0]

  specVec = []
  for i, mapType in enumerate(mapTypes):
    specVec.append(str(specArr[i][0]))
    specVec.append(str(specArr[i][1]))

  specStr = ','.join(specVec)  
  return specVec      
  
