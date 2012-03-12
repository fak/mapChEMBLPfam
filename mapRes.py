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
  
  
def findPrimers(pdbDict, pfamDict):
  import numpy as np
  
  primers = {}
  for target in pdbDict.keys():
    for cmpdId in pdbDict[target].keys():
      tmpDict = {}

      domains = {}
      for domain in pfamDict[target]['domains']:
        domains[domain] = 0
      for domain in domains:
        ndom = len([x for x in pdbDict[target][cmpdId]['domain'] if x == domain])
        nall = len(pdbDict[target][cmpdId]['domain'])
        ratio = np.true_divide(ndom, nall)
        if ndom > 4 and ratio >= 0.3:
          tmpDict[domain] = ratio
      # print i, len(pdbDict[target][cmpdId]['domain']), len(pdbDict[target][cmpdId]['pdb']), target, cmpdId
      primer= ' $$$ '.join(tmpDict.keys())
      ratios = tmpDict.values()
      pdbs = {}
      for pdb in pdbDict[target][cmpdId]['pdb']:
        pdbs[pdb] = 0

      try:
        primers[primer]['ratios'].append(ratios)
        primers[primer]['pdb'].append(pdbs.keys())
        primers[primer]['cmpdId'].append(cmpdId)
      except KeyError:
        primers[primer] = {}
        primers[primer]['pdb'] = []
        primers[primer]['cmpdId'] = []
        primers[primer]['ratios'] = []
        primers[primer]['cmpdId'].append(cmpdId)
        primers[primer]['ratios'].append(ratios)
        primers[primer]['pdb'].append(pdbs.keys())
      
  return primers            
