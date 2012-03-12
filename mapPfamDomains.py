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

  ## Load the pdbDict.
  infile = open('data/pdbDict_chembl13.pkl','r')
  pdbDict = pickle.load(infile)
  infile.close()

  ###    ###    ###     ###  
  ###    ###    ###     ###
  ###    ###    ###     ###
  
  ### Get primers from pdbe.
  import mapRes
  pdbDict = mapRes.mapRes(pdbDict, pfamDict, release)
  primers = prime.findPrimers(pdbDict, pfamDict)
  hierDict = prime.hier(primers)
  
  # Assign protein targets to the primers.
  primers = prime.mapTargets(chemblTargets, primers, hierDict)  
  
  # Find conflicts.
  conflicts = findConflicts.findConflicts(primers)
  confTargets = findConflicts.confTargets(conflicts)
  
  #Populate the primer classes with active compounds. Exempt compounds occuring in conflicting cases.
  mpdct = mpdct.mpdct(confTargets, threshold, release, user, pword, host, port)
   
  ## Construct the propDict for targets with one domain and manually 
  ## insert/delete ligands.
  import feedPropDict
  import parse
  blacklist = parse.col2list('data/removeLigands.txt',1, False)  
  propDict = {}
  propDict = feedPropDict.dictionary(mpdct, propDict, blacklist, 'first')
  

  ## Generate output files for Sam.
  #import toSam
  #toSam.toSam(conflicts, threshold, user, pword, host, release, port)
  ## Use input files from Sam to map the conflicts.
  #import fromSam
  #conf = fromSam.fromSam(conflicts, threshold, user, pword, host, release, port)

  ## Insert data for resolved conflicts.
  #import feedPropDict
  #propDict = feedPropDict.dictionary(conf, propDict, blacklist,'conflict')

  ## Export the mapping to a mySQL table.
  import export
  import pickle
  outfile = open('data/propDict_%s.pkl' %release, 'w')
  pickle.dump(propDict, outfile)
  export.exportMapsMySQL(propDict, release, user, pword, host, port)
  export.exportConflsMySQL(conflicts, release ,user, pword, host, port)
  
