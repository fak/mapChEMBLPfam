"""
  Function:  generateMolDict
  --------------------
  Retrieves a list of PDBe identifiers and matches them to molregnos.
  
  momo.sander@ebi.ac.uk
"""                                       
def generateMolDict(release, user, pword, host, port):

  import queryDevice

  cmpds = queryDevice.queryDevice("SELECT molregno, compound_key \
  FROM compound_records WHERE src_id = 6", release, user, pword, host, port)

  molDict = {}
  for cmpd in cmpds:
    molDict[cmpd[0]] = cmpd[1]

  return molDict

