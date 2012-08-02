"""
Function:  uniprot(bsDict, pfamDict release) 

Checks whether the binding site is outside all domains and adds a boolean to 
pfamDict[target]['prediction'] for each listed residue

momo.sander@googlemail.com
"""

def uniprot(bsDict, pfamDict, release):

  for target in bsDict.keys():
    try:
      pfamDict[target]
    except KeyError:
      #print 'No prediction for %s'%target
      del bsDict[target]
      continue
    bsDict[target]['within'] = []
    for pos in bsDict[target]['positions']:
      i = 0
      pred = False
      pos = int(pos)
      for domain in pfamDict[target]['domains']:
        
        start = pfamDict[target]['start'][i]
        end = pfamDict[target]['end'][i]
        i +=1
        #print pos, start, end
        if pos >= start and pos <= end:
          pred = True

  
      bsDict[target]['within'].append(pred)
      
  return bsDict


"""
  Function:  pdbe
  --------------------
  Check if positions are within any given domain. 
  
  momo.sander@ebi.ac.uk
"""                        
def pdbe(pdbDict,pfamDict, release):

  for target in pfamDict.keys():
    try:
      pdbDict[target]
    except KeyError:
      continue
    for cmpdId in pdbDict[target].keys():
      for pos in pdbDict[target][cmpdId]['position']:
        pos = int(pos)
        i = 0
        pred = False
        for domain in pfamDict[target]['domains']:
          
          start =pfamDict[target]['start'][i]
          end = pfamDict[target]['end'][i]
          i+=1
          #print 'Interacting residue, start, end:', pos, start, end, '\n' 
          if pos >= start  and pos <= end :
            pred = True
        try:
          pdbDict[target][cmpdId]['within'].append(pred)
        except KeyError:
          pdbDict[target][cmpdId]['within'] = []
          pdbDict[target][cmpdId]['within'].append(pred)

  return pdbDict

"""
  Function:  pdbePredicted
  --------------------
  Check if positions are within the predicted domain. 
  
  momo.sander@ebi.ac.uk
"""  

 
def pdbePredicted(pdbDict, intacts, uniDict):
   
  for intact in intacts:
    target = intact[0]
    start = intact[3]
    end = intact[4]
    mapType = intact[5]
    molregno = intact[2]
    chemblId = intact[6]
    try:
      code = uniDict[chemblId][0]
    except KeyError:  
      continue
    try:
      pdbDict[target][code]
    except KeyError:
      continue

    pdbDict[target][code]['maptype'] = mapType
    for pos in pdbDict[target][code]['position']:
      pos = int(pos)
      i = 0
      pred = False
      if pos >= start  and pos <= end :
        pred = True
      try:
        pdbDict[target][code]['prediction'].append(pred)
      except KeyError:
        pdbDict[target][code]['prediction'] = []
        pdbDict[target][code]['prediction'].append(pred)
    

  return pdbDict 



"""
  Function:  uniprotPredicted
  --------------------
  Check if positions are within the predicted domain. 
  
  momo.sander@ebi.ac.uk
"""


def uniprotPredicted(bsDict, intacts):

  for intact in intacts:
    target = intact[0]
    start = intact[3]
    end = intact[4]
    mapType = intact[5]
    
    try:
      bsDict[target]
    except KeyError:
      continue

    bsDict[target]['maptype'] = mapType
    for pos in bsDict[target]['positions']:
      pos = int(pos)
      i = 0
      pred = False
      if pos >= start  and pos <= end :
        pred = True
      try:
        bsDict[target]['prediction'].append(pred)
      except KeyError:
        bsDict[target]['prediction'] = []
        bsDict[target]['prediction'].append(pred)


  return bsDict





