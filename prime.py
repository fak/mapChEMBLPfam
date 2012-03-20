"""
Function:  prime
  --------------------
  Collection of functions dealing with the primers dictionary.
  
  momo.sander@ebi.ac.uk
"""                        
def findPrimers(pdbDict, pfamDict):
  import numpy as np
  
  primers = {}
  for target in pdbDict.keys():
    if target not in pfamDict.keys():
      continue
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


"""
  Function:  hier
  --------------------
  # Create a look-up to deal with primer hierarchy.
  
  momo.sander@ebi.ac.uk
""" 
def hier(primers):

  
  lkp = {}
  for primer1 in sorted(primers.keys()):
    elements1 = primer1.split(' $$$ ')
    for primer2 in sorted(primers.keys()):
      elements2 = primer2.split(' $$$ ') 
      count = 0
      for element in elements1:
        if element in elements2:
          count +=1
      if count >= len(elements1) and primer1 != primer2:
        try:  
          lkp[primer1].append(primer2)
        except KeyError:
          lkp[primer1] = [primer2]
  
  return lkp
  

"""
  Function:  mapTargets
  --------------------
  # Assign protein targets to the primers
  
  momo.sander@ebi.ac.uk
"""
def mapTargets(chemblTargets, pfamDict,  primers, hierDict):
  for primer in primers.keys():
    primers[primer]['targets'] = {}
  for target in chemblTargets:
    try:
      pfamDict[target]
    except KeyError:
      continue
      
    for primer in primers.keys():
      alt = False
      try:
        altPrimers = hierDict[primer]
      except KeyError:
        altPrimers = []
        
      for altprim in altPrimers:
        altDomains = altprim.split(' $$$ ')
        lkp = {}
        for domain in altDomains:
          if domain in pfamDict[target]['domains']:
            lkp[domain] = 0
        if sorted(altDomains) == sorted(lkp.keys()):
          alt = True
          #print primer, hierDict[primer]
          
      if not alt:
        domains = primer.split(' $$$ ')
        lkp = {}  
        for domain in domains:
          if domain in pfamDict[target]['domains']:
            lkp[domain] = 0
        if sorted(domains) == sorted(lkp.keys()):
          primers[primer]['targets'][target]=0
            
  return primers
