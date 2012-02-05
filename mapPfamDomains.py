"""
  Function:  mapPDsa
  --------------------
  Carry out the mapping and save results.
  
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
  inFile = open('data/protCodPfamDict_%s.pkl' %release, 'r')
  pfamDict = pickle.load(inFile)
  inFile.close()    

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
  ## use Sam's procedure.

  ## Insert data for conflicts and multi domain proteins
  import feedPropDict
  propDict = feedPropDict.dictionary(multi, propDict, blacklist, 'multi')
  #propDict = feedPropDict.dictionary(confDict, propDict, blcklist,'conflict')

  ## Export the mapping to a mySQL table.
  import export
  import pickle
  outfile = open('data/propDict_%s' %release, 'w')
  pickle.dump(propDict, outfile)
  export.exportMapsMySQL(propDict, release, user, pword, host, port)
  export.exportConflsMySQL(conflicts, release ,user, pword, host, port)
  
