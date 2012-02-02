
"""
  Function:  filterForTarget
  --------------------
"""

def filterForTarget(ligands, threshold):

  print threshold
  import pickle
  import math
  import numpy as np
  print 'filtering ligands at threshold:', threshold
  ligandList = []
  for data in ligands:   
    pAfnty=0  
    molregno = data[0]
    standardValue = data[1]
    
    try:
      standardValue = float(standardValue)
    except:
      continue
      
    standardValue = float(standardValue) 
    if standardValue <=0:  
      continue
    standardType = data[2]
    standardUnit = data[3]
    smiles = data[4]
    relation = data[5]
    docId = data[6]
    if standardType in ['Ki','IC50', 'EC50'] and standardUnit == 'nM' and relation == '=':
      pAfnty = (-1) * math.log10(standardValue/float(1000000000))
    elif standardType in ['-Log Ki','pKd', 'pA2', 'pI', 'pKi']:
      pAfnty = standardValue
    
    else:
      continue
   
    if pAfnty <= threshold:
      continue
    
    else:
      ligandList.append([smiles,pAfnty,molregno, docId])
  return ligandList

 
