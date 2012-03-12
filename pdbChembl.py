"""
  Function:  query
  --------------------
  Retrieves all binding site residues for interactions in ChEMBL.
  
  momo.sander@ebi.ac.uk
"""                                       
                                     
def query(release, user, pword, host, port): 
                           

  ## Load the uniDict.
  import parseUniChem
  uniDict = parseUniChem.parse('data/unichemMappings.txt')

  ## Create the intactDict.          
  import getIntactDict
  intactDict =  getIntactDict.getIntacts(uniDict, release, user, pword, host, port)
  
  ## Load the pdbDict.
  import queryPDB
  pdbDict = queryPDB.queryPDB(uniDict, intactDict, release)
  
  return pdbDict                   
