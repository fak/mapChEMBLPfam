"""
Function:  parseUniChem
  --------------------
  Carry out the mapping and save results.

  Author:
  Felix Kruger
  momo.sander@googlemail.com
"""                              
def parse(infile): 
  
  infile = open(infile, 'r')
  lines = infile.readlines()  
  uniDict = {}
  
  for line in lines[1:]:
      if len(line) == 1:
          continue
      elements  = line.split("'--->'")
      pdb = elements[0].lstrip("'")      
      chembl = elements[1].rstrip("'\n")
      if len(chembl) ==0 or len(pdb) == 0:
          continue
      try:
        uniDict[chembl].append(pdb)
      except KeyError:
        uniDict[chembl] = []
        uniDict[chembl].append(pdb)
  
  return uniDict   
