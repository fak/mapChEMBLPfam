"""
  Function:  singleDomains
  --------------------
  Identify targets with one binding site containing domain and at least one other domain and train a generic model
  momo.sander@ebi.ac.uk
"""                              

def singleDomains(pfamDict, chemblTargets, release):
  import pickle
  import getLigands
  import filterForTarget
 
  single = {}
  for target in chemblTargets:
    if target in pfamDict:
      if len(pfamDict[target]['domains']) == 1:
        domain = pfamDict[target]['domains'][0]
        ligands = getLigands.getLigandsForTarget(target, release)
        ligands = filterForTarget.filterForTarget(ligands, threshold)
        for ligand in ligands:
          smiles = ligand[0]
          aff = ligand[1]
          molregno = ligand[2]
          actId = ligand[3]
          try:
            single[domain][molregno]['pAfnty'].append(aff)
            single[domain][molregno]['target'].append(target)
	        single[domain][molregno]['actId'].append(actId)
            single[domain][molregno]['smiles']=smiles
          except KeyError:
            try:
              single[domain][molregno] = {}
              single[domain][molregno]['pAfnty']=[]
              single[domain][molregno]['target']=[]
	          single[domain][molregno]['actId'] = []
              single[domain][molregno]['smiles']=smiles
              single[domain][molregno]['pAfnty'].append(aff)
              single[domain][molregno]['target'].append(target)
              single[domain][molregno]['actId'].append(actId)
            except KeyError:
              single[domain]={}
              single[domain][molregno] = {}
              single[domain][molregno]['pAfnty']=[]
              single[domain][molregno]['target']=[]
              single[domain][molregno]['actId'] = []
              single[domain][molregno]['smiles']=smiles
              single[domain][molregno]['pAfnty'].append(aff)
              single[domain][molregno]['target'].append(target)
              single[domain][molregno]['actId'].append(actId)
  
  outfile = open('data/singleDict_pKi%s_%s.pkl' %(int(threshold), release),'w')
  pickle.dump(single, outfile)
  outfile.close()  

  return single


                                          
                                                                         
    
