"""
  Function:  mapPDs
  --------------------
  Carry out the mapping and save results.

  Author:
  Felix Kruger
  momo.sander@googlemail.com
"""                              
def mapPDs(release, user, pword, host, port): 

  ## Set the threshold.
  import numpy as np
  threshold = -np.log10(50*10**(-6))
  
  ## Get a list of all ChEMBL targets.
  import getUniprotTargets
  chemblTargets = getUniprotTargets.getUniprotTargets(release, user, pword, host, port)

  ## Load the pfamDict.
  import pickle
  infile = open('data/protCodPfamDict_%s.pkl' %release, 'r')
  pfamDict = pickle.load(inFile)
  infile.close()    

  ## Load the uniDict.
  import parseUniChem
  uniDict = parseUniChem.parse('data/unichemMappings.txt')

  ## Load the intactDict.          
  import getIntactDict
  intactDict =  getIntactDict.getIntacts(uniDict, release, user, pword, host, port)
  
  ## Load the pdbDict.
  import queryPDB
  pdbDict = queryPDB.queryPDB(uniDict, intactDict, release)

  infile = open('data/pdbDict_chembl13.pkl','r')
  pdbDict = pickle.load(infile)
  infile.close()

  ###    ###    ###     ###  
  ###    ###    ###     ###
  ###    ###    ###     ###
  
  ### Get primers from pdbe.
  import mapRes
  pdbDict = mapRes.mapRes(pdbDict, pfamDict, release)
  primers = mapRes.findPrimers(pdbDict, pfamDict)

  # Create a look-up to deal with primer hierarchy.
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
  
  hierDict = lkp

  
  # Find all primers that are contained within another primer and clear shared targets.      
  #p primers t,primer1 in enumerate(sorted(primers.keys())[:-1]):
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

  # Find conflicts.
  conflicts = {}
  for i,primer1 in enumerate(sorted(primers.keys())[:-1]):
    for primer2 in sorted(primers.keys())[i+1:]:
      if primer1 == primer2: print i, primer1
      targets1 = primers[primer1]['targets'].keys()
      for target in targets1:
        if target in primers[primer2]['targets'].keys():
          conflict_string = '%%%'.join(sorted([primer1, primer2]))
          try:
            conflicts[conflict_string].append(target)
          except KeyError:
            conflicts[conflict_string] = []
            conflicts[conflict_string].append(target)


  #Populate the primer classes with active compounds. Exempt compounds occuring in conflicting cases.
  

sys.exit()



  ## Get ligands for targets with single domains.
  import singleDomain 
  single = singleDomain.singleDomains(pfamDict, chemblTargets, threshold, release, user, pword, host, port)

  ## Construct the propDict for targets with one domain and manually 
  ## insert/delete ligands. From this, we also derive the list of winners, ie. 
  ## domains binding ligands with at least micromolar affinity. 
  import feedPropDict
  import parse
  blacklist = parse.col2list('data/removeLigands.txt',1, False)  
  propDict = {}
  propDict = feedPropDict.dictionary(single, propDict, blacklist, 'single')
  propDict = feedPropDict.addLigs(propDict,'manual') 
  
  ## Extract a list of validated domains.
  valid = propDict.keys() 

  ## Identify targets with one binding site containing domain and at least one
  ## other domain.
  import multiDomain
  multi = multiDomain.multiDomain(pfamDict, chemblTargets, valid, threshold, release, user, pword, host, port)

  ## Deal with targets that have more than one binding site containing 
  ## domain. Interface with Sam's text mining process.
  import findConflicts
  conflicts = findConflicts.findConflicts(pfamDict, valid, chemblTargets)
  ## Generate output files for Sam.
  import toSam
  #toSam.toSam(conflicts, threshold, user, pword, host, release, port)
  ## Use input files from Sam to map the conflicts.
  import fromSam
  conf = fromSam.fromSam(conflicts, threshold, user, pword, host, release, port)

  ## Insert data for conflicts and multi domain proteins
  import feedPropDict
  propDict = feedPropDict.dictionary(multi, propDict, blacklist, 'multi')
  propDict = feedPropDict.dictionary(conf, propDict, blacklist,'conflict')

  ## Export the mapping to a mySQL table.
  import export
  import pickle
  outfile = open('data/propDict_%s.pkl' %release, 'w')
  pickle.dump(propDict, outfile)
  export.exportMapsMySQL(propDict, release, user, pword, host, port)
  export.exportConflsMySQL(conflicts, release ,user, pword, host, port)
  
