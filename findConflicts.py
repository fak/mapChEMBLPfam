"""
  Function:  findConflicts
  --------------------
  searches for proteins containing more than one ligand binding domain and 
  identifies the combinations of domains that are in conflict with each other... 
  
  momo.sander@ebi.ac.uk
"""                                       
def findConflicts(primers): 

  conflicts = {}
  for i,primer1 in enumerate(sorted(primers.keys())[:-1]):
    for primer2 in sorted(primers.keys())[i+1:]:
      if primer1 == primer2: print i, primer1
      targets1 = primers[primer1]['targets'].keys()
      for target in targets1:
        if target in primers[primer2]['targets'].keys():
          conflict_string = '%%%'.join(sorted([primer1, primer2]))
          try:
            conflicts[conflict_string].append(target)
          except KeyError:
            conflicts[conflict_string] = []
            conflicts[conflict_string].append(target)

  return conflicts
  
  
def confTargets(conflicts): 
  confTargets = []
  for conf in conflict.keys():
    for target in conflicts[conf]:
      confTargets.append(target)

  return confTargets
