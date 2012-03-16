"""
  Function:  findConflicts
  --------------------
  searches for proteins containing more than one ligand binding domain and 
  identifies the combinations of domains that are in conflict with each other... 
  
  momo.sander@ebi.ac.uk
"""                                       
def findConflicts(pfamDict, winners, chemblTargets): 

  conflicts = {}

  for target in chemblTargets:
    candidates = {}
    if target in pfamDict:
      for domain in pfamDict[target]['domains']:
        domain = str(domain)      
        if domain in winners:
          candidates[domain]=0
        #else:
         # print domain, ' is not a winner'
    else:
      continue         
    if len(candidates.keys()) > 1:
      conflict_string = '%%%'.join(sorted(candidates.keys()))
      try: 
        conflicts[conflict_string].append(target)
      except KeyError:
        conflicts[conflict_string] = []
        conflicts[conflict_string].append(target)

  return conflicts
