"""
Function:  uniprotChembl 


momo.sander@googlemail.com
"""

def query(release, user, pword, host, port):
  
  import getBindingSitesUniprot
  import getUniprotTargets


  ## Get all protein target s from ChEBML.                                       
  chemblTargets = getUniprotTargets.getUniprotTargets(release, user, pword, host, port)

  ## Get binding sites for each target.
  bsDict = getBindingSitesUniprot.getBindingSites(chemblTargets, release)
  print 'number of targets with binding site information', len(bsDict.keys())

