"""
  Function:  pdbChembl
  --------------------
  Retrieves all binding site residues for interactions in ChEMBL.
  
  fkrueger@ebi.ac.uk
"""                                       
                                     
def query(release, user, pword, host, port): 
                           

  import getIntactDict
  import queryPDB
  import coordMap
  import parseUniChem
  ## Create the pdbDict which stores all the relevant information. The moldict
  ## is a look-up molDict[molregno] = cmpdId. Coordmap is a dictionary with the
  ## coordinate mappings between pdbe and uniprot.

  uniDict = parseUniChem.parse('data/unichemMappings.txt')
  coordMap = coordMap.coordMap()
  intactDict =  getIntactDict.getIntacts(uniDict, release, user, pword, host, port)
  pdbDict = queryPDB.queryPDB(uniDict, intactDict, coordMap,  release)
  return pdbDict                   
