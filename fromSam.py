"""
Function: fromSam
  --------------------
  momo.sander@googlemail.com
"""                    

def fromSam(conflicts, threshold, user, pword, host, release, port):
  
  import parse
  import getLigands
  import filterForTarget
  import queryDevice 

  conf = {}
  for conflict in conflicts.keys():
    confStr = conflict.split('%%%')
    confStr = '_'.join(confStr)
    try:
      predMaps = parse.parse2col( 'data/%s.processed'%confStr, True, 3,4) 
    except IOError:
      continue

    for act in predMaps.keys():
      if predMaps[act] == 'None':
        continue
      else:
        domain = predMaps[act]       
        ligands = getLigands.getLigandsForActivity(act, release, user, pword, host, port)
        target = queryDevice.queryDevice("SELECT protein_accession FROM target_dictionary td JOIN assay2target a2t ON a2t.tid = td.tid JOIN activities act ON a2t.assay_id = act.assay_id WHERE activity_id = %s"%act, release, user, pword, host, port) 
        
        ligands = filterForTarget.filterForTarget(ligands, threshold)
        for ligand in ligands:
          smiles = ligand[0]
          aff = ligand[1]
          molregno = ligand[2]
          actId = ligand[3]
          try:
            conf[domain][molregno]['pAfnty'].append(aff)
            conf[domain][molregno]['target'].append(target)
            conf[domain][molregno]['actId'].append(actId)
            conf[domain][molregno]['smiles']=smiles
          except KeyError:
            try:
              conf[domain][molregno] = {}
              conf[domain][molregno]['pAfnty']=[]
              conf[domain][molregno]['target']=[]
              conf[domain][molregno]['actId'] = []
              conf[domain][molregno]['smiles']=smiles
              conf[domain][molregno]['pAfnty'].append(aff)
              conf[domain][molregno]['target'].append(target)
              conf[domain][molregno]['actId'].append(actId)
            except KeyError:
              conf[domain]={}
              conf[domain][molregno] = {}
              conf[domain][molregno]['pAfnty']=[]
              conf[domain][molregno]['target']=[]
              conf[domain][molregno]['actId'] = []
              conf[domain][molregno]['smiles']=smiles
              conf[domain][molregno]['pAfnty'].append(aff)
              conf[domain][molregno]['target'].append(target)
              conf[domain][molregno]['actId'].append(actId)
  return conf   

