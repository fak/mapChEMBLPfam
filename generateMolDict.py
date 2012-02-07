"""
  Function:  generateMolDict
  --------------------
  Creates a mapping between Molregnos and PDBe identifiers. 
  
  momo.sander@ebi.ac.uk
"""                                       
def gen(intacts, release):
  import queryDevice
  
  molDict = {}
  for tup in intacts:
    molDict[tup[2]] = 0

  for molregno in molDict.keys():
    codes = queryDevice.queryDevice("SELECT compound_key FROM compound_records WHERE molregno = %s AND src_id =6" % molregno, release) 
    for code in codes:
      print code[0]
      molDict[molregno] = code[0]

  return molDict


