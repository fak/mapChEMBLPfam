"""
Function:  redThread 


momo.sander@googlemail.com
"""

def redThread(release):
  
  import getBindingSitesUniprot
  import getUniprotTargets


  ## Get all protein target s from ChEBML.                                       
  chemblTargets = getUniprotTargets.getUniprotTargets(release)

  ## Get binding sites for each target.
  bsDict = getBindingSitesUniprot.getBindingSites(chemblTargets, release)
  print 'number of targets with binding site information', len(bsDict.keys())

