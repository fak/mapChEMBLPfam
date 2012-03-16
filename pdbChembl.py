"""
  Function:  query
  --------------------
  Retrieves all binding site residues for interactions in ChEMBL.
  
  momo.sander@ebi.ac.uk
"""                                       
                                     
def query(release, user, pword, host, port): 
                           

  import makeIntactDict
  import queryPDB
  import coordMap
 
  ## Create the pdbDict which stores all the relevant information. The moldict
  ## is a look-up molDict[molregno] = cmpdId. Coordmap is a dictionary with the
  ## coordinate mappings between pdbe and uniprot.

    
  coordMap = coordMap.coordMap()
  intactDict =  getIntactDict.getIntacts(uniDict, release, user, pword, host, port)
  pdbDict = queryPDB.queryPDB(uniDict, intactDict, coordMap,  release)
  return pdbDict                   
