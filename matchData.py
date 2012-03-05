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
    bsDict[target]['prediction'] = []
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

  
      bsDict[target]['prediction'].append(pred)
      
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
          pdbDict[target][cmpdId]['prediction'].append(pred)
        except KeyError:
          pdbDict[target][cmpdId]['prediction'] = []
          pdbDict[target][cmpdId]['prediction'].append(pred)

  return pdbDict

"""
  Function:  pdbePredicted
  --------------------
  Check if positions are within the predicted domain. 
  
  momo.sander@ebi.ac.uk
"""  

 
def pdbePredicted(pdbDict, intacts, molDict, release, mapType):
  import queryDevice

  predList = [] 
  for intact in intacts:
    target = intact[0]
    start = intact[3]
    end = intact[4]
    molregno = intact[2]
    try:
      code = molDict[molregno]
    except KeyError:  
      continue
    try:
      pdbDict[target][code]
    except KeyError:
      continue
      #for dummyTarget in pdbDict.keys():
       # if code in pdbDict[dummyTarget].keys():
        #  target = dummyTarget
        #  continue
         
    preds = ['%s_%s_%s'%(target, molregno, code)]
    for pos in pdbDict[target][code]['position']:
      pos = int(pos)
      i = 0
      pred = False
      #print 'Interacting residue, start, end:', pos, start, end, '\n' 
      if pos >= start  and pos <= end :
        pred = True
      preds.append(pred)
    predList.append(preds)

  return predList  


