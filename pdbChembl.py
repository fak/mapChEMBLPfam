"""
  Function:  query
  --------------------
  Retrieves all binding site residues for interactions in ChEMBL.
  
  momo.sander@ebi.ac.uk
"""                                       
                                     
def query(release, user, pword, host, port): 
                           

  import makeIntactDict
  import queryPDB

  ## Create the pdbDict which stores all the relevant information. The moldict
  ## is a look-up molDict[molregno] = cmpdId.
  
  (intactDict, molDict) = makeIntactDict.mkIntactDictAllPDBs(release,user,pword, host, port)
  pdbDict = queryPDB.queryPDB(molDict, intactDict, release)
  
  return pdbDict                   
