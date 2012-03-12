"""
  Function:  mapRes 
  --------------------
  Check if positions are within any given domain. 
  
  momo.sander@ebi.ac.uk
"""                        
def mapRes(pdbDict,pfamDict, release):

  for target in pfamDict.keys():
    try:
      pdbDict[target]
    except KeyError:
      continue
    for cmpdId in pdbDict[target].keys():
      pdbDict[target][cmpdId]['domain']=[]
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
            pred = domain
            pdbDict[target][cmpdId]['domain'].append(pred)
            break
        if not pred:
          pdbDict[target][cmpdId]['domain'].append(pred)
            
  return pdbDict
