"""
  Function:  prime
  --------------------
  Collection of functions dealing with the primers dictionary.
  
  momo.sander@ebi.ac.uk
"""                       
def exportPrimers(primers):
  outF = open('data/markdownFullCombis_%s.md' % release ,'w')
  outF.write('|domain combination|pdb accession|ratios|\n')
  outF.write('|:-----------|------------:|:------------:|\n')
  
  out = open('data/markdownCombis_%s.md' % release ,'w')
  out.write('|domain combination|pdb accession|# ChEMBL targets|\n')
  out.write('|:-----------|------------:|:------------:|\n')
  for primer in sorted(primers.keys()):
    els = primer.split(" $$$ ")
    if len(els) < 2:
      continue
    pString = ", ".join(els)
    pdbs = {}
    for i, pdb in enumerate(primers[primer]['pdb']):
      pdbs[pdb[0]] = 0
      ratios = primers[primer]['ratios'][i]
      rstring = ', '.join(['%.2f' % x for x in ratios])
      outF.write("|%s|%s|%s|\n"%(pString, pdb[0], rstring))
  
    pdbString = ", ".join(pdbs.keys()[:min(3,len(pdbs.keys()))])
    ntargs = len(primers[primer]['targets'])
    out.write("|%s|%s|%s|\n"%(pString, pdbString, ntargs))
  out.close()
  outF.close()

 
def longPfamDict(pfamDict):
  longPD = {}
  for accession in pfamDict.keys():
    longPD[accession] = {}
    for i, domain in enumerate(pfamDict[accession]['domains']):
      start = pfamDict[accession]['start'][i]
      end = pfamDict[accession]['end'][i]
      for j in range(start, end):
        longPD[accession][j] = domain
  return longPD


def findPrimers(pdbDict, longPD, min_res, min_ratio):
  import numpy as np
  primers = {}
  for target in pdbDict.keys():
    domains = {}
    for domain in longPD[target].values():
      try:
        domains[domain] += 1
      except KeyError:
        domains[domain] = 1
    for cmpdId in pdbDict[target].keys():
      tmpDict = {}
      for domain in domains:
        ndom = len([x for x in longPD[target].keys() if longPD[target][x] == domain and x in pdbDict[target][cmpdId]['position']])
        nall = len(pdbDict[target][cmpdId]['position'])
        ratio = np.true_divide(ndom, nall)
        if ndom >= min_res and ratio >= min_ratio:
          tmpDict[domain] = ratio
      if len(tmpDict.keys()) <2:
        continue
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
        primers[primer]['pdb'] = [pdbs.keys()]
        primers[primer]['cmpdId'] = [cmpdId]
        primers[primer]['ratios'] = [ratios]
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
def mapTargets(chemblTargets, primers, hierDict):
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
          print primer, hierDict[primer]
          
      if not alt:
        domains = primer.split(' $$$ ')
        lkp = {}  
        for domain in domains:
          if domain in pfamDict[target]['domains']:
            lkp[domain] = 0
        if sorted(domains) == sorted(lkp.keys()):
          try:
            primers[primer]['targets'][target]=0
          except KeyError: 
            primers[primer]['targets']={}
            primers[primer]['targets'][target]=0
            
  return primers
