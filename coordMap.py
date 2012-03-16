"""
  Function:  coordMap
  --------------------
  Parse the sift mapping and creates a dictionary with the coordinate mappings.
  
  momo.sander@ebi.ac.uk
"""

def coordMap():

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

  return coordMap
