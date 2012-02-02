"""
 Generate a list of all targets that are to be fed into the getPfamDomain 
 procedure.
  --------------------
  Author:
  Felix Kruger
  momo.sander@googlemail.com
"""        



def getAllTargets(humanTargets, chemblTargets):
  tDict = {}

  for target in chemblTargets:
    tDict[target] = 'chembl'
  for target in humanTargets:
    tDict[target] = 'human'

  return tDict
