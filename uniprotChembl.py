"""
Function:  uniprotChembl 


momo.sander@googlemail.com
"""

def query(release, user, pword, host, port):
  
  import queryUniprot
  import getUniprotTargets


  ## Get all protein targets from ChEBML.                                       
  chemblTargets = getUniprotTargets.getUniprotTargets(release, user, pword, host, port)

  ## Get Uniprot binding site annotation for each target.
  uniDict = queryUniprot.getBindingSites(chemblTargets, release)
  print 'number of targets with binding site information', len(uniDict.keys())

