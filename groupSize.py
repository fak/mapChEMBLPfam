"""
  Function:  groupSize
  --------------------
  calculates the group sizes for each mapping case 
  
  momo.sander@ebi.ac.uk
"""                              
def groupSize(chemblTargets, pfamDict): 
  import numpy as np

  cMult = []
  cSing = []
  cNone = []
  cConf = []
  
  winners = []
  for target in chemblTargets:
    try:
      count = pfamDict[target]['countUnique']
      if count == 1:
        winners.append(pfamDict[target]['domains'][0])
    except KeyError:
      print 'not in Pfam: ', target

  for target in chemblTargets:
    try:
      domains = pfamDict[target]['domains']
      if len(domains) > 1:
        wins = {}
        for domain in domains:
          if domain in winners:
            wins[domain] = 0
        if len(wins.keys())>1:
          cConf.append(target)
        elif len(wins.keys()) == 1:
          cMult.append(target)
        elif len(wins.keys()) == 0:
          cNone.append(target)
      elif len(domains) == 1:
        cSing.append(target)
    except KeyError:
      print 'not in Pfam: ', target

  groups = [str(len(cSing)), str(len(cNone)), str(len(cMult)), str(len(cConf))]
  args  = ','.join(groups)
  return args
  

def groupSizeMap(chemblTargets, release, user, pword, host, port): 
  import queryDevice

  multi = queryDevice.queryDevice("SELECT DISTINCT protein_accession FROM map_pfam WHERE mapType = 'multi'", release, user, pword, host, port)
  
  single = queryDevice.queryDevice("SELECT DISTINCT protein_accession FROM map_pfam WHERE mapType = 'single'", release, user, pword, host, port)
  
  conflict = queryDevice.queryDevice("SELECT DISTINCT protein_accession FROM map_pfam WHERE mapType = 'conflict'", release, user, pword, host, port)

 
  return single,multi, conflict
 

def actSizeMap(chemblTargets, release, user, pword, host, port):
  import queryDevice

  multi = queryDevice.queryDevice("SELECT DISTINCT activity_id FROM map_pfam WHERE mapType = 'multi'", release, user, pword, host, port)

  single = queryDevice.queryDevice("SELECT DISTINCT activity_id FROM map_pfam WHERE mapType = 'single'", release, user, pword, host, port)

  conflict = queryDevice.queryDevice("SELECT DISTINCT activity_id FROM map_pfam WHERE mapType = 'conflict'", release, user, pword, host, port)


  return  single,multi, conflict
 
