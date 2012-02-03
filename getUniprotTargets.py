"""
Function:  getUniprotTargets

  gets all protein accession with db_source swiss-prot or trembl from the Chembl 
  target dictionary
  
  momo.sander@googlemail.com
"""

def getUniprotTargets(release, user, pword, host, port):

  import queryDevice
  
  rawtargets = queryDevice.queryDevice("SELECT protein_accession\
    FROM target_dictionary WHERE db_source IN('SWISS-PROT', 'TREMBL')", release, user, pword, host, port)

  targets= []
  tids = []
  for target in rawtargets:
    targets.append(target[0])


  

  return targets

