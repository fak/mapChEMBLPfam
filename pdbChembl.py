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

  #need to insert a script that parses the sift mapping and creates a dictionary with the coordinate mappings.
  infile = open('data/pdb_chain_uniprot.csv', 'r')
  lines = infile.readlines()
  coordMap = {}
  for line in lines[1:]:
    elements = line.split(',')
    pdb = elements[0]
    chain = elements[1]
    uniprot = elements[2]
    try:
      pdbStart = int(elements[5])
      uniprotStart = int(elements[7])
    except ValueError:
      continue
    offset = uniprotStart - pdbStart
    try: 
      coordMap[pdb][chain] = offset
    except KeyError:
      coordMap[pdb] = {}
      coordMap[pdb][chain] = offset 
    
  intactDict =  getIntactDict.getIntacts(uniDict, release, user, pword, host, port)
  pdbDict = queryPDB.queryPDB(uniDict, intactDict, coordMap,  release)
  return pdbDict                   
