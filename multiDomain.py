"""
  Function:  multiDomain
  --------------------
  Identify targets with one binding site containing domain and at least one other domain and train a generic model
  momo.sander@ebi.ac.uk
"""                              
def multiDomain(pfamDict, chemblTargets, winners, threshold, release, user, pword, host, port):

  import getLigands
  import filterForTarget
  import pickle

 
  multi = {}
  for target in chemblTargets:
    candidates = {}
    if target in pfamDict:
      for domain in pfamDict[target]['domains']:
        domain = str(domain)      
        if domain in winners:
          candidates[domain]=0
        
    if len(candidates.keys()) == 1 and len(pfamDict[target]['domains']) >1:
      ligands = getLigands.getLigandsForTarget(target, release, user, pword, host, port)
      ligands = filterForTarget.filterForTarget(ligands, threshold)
      domain = candidates.keys()[0]
      for ligand in ligands:
        smiles = ligand[0]
        aff = ligand[1]
        molregno = ligand[2]
	actId = ligand[3]
        try:
          multi[domain][molregno]['pAfnty'].append(aff)
          multi[domain][molregno]['target'].append(target)
	  multi[domain][molregno]['actId'].append(actId)
          multi[domain][molregno]['smiles']=smiles
        except KeyError:
          try:
            multi[domain][molregno] = {}
            multi[domain][molregno]['pAfnty']=[]
            multi[domain][molregno]['target']=[]
            multi[domain][molregno]['actId'] = []
            multi[domain][molregno]['smiles']=smiles
            multi[domain][molregno]['pAfnty'].append(aff)
            multi[domain][molregno]['target'].append(target)
	    multi[domain][molregno]['actId'].append(actId)
          except KeyError:
            multi[domain]={}
            multi[domain][molregno] = {}
            multi[domain][molregno]['pAfnty']=[]
            multi[domain][molregno]['target']=[]
	    multi[domain][molregno]['actId'] = []
            multi[domain][molregno]['smiles']=smiles
            multi[domain][molregno]['pAfnty'].append(aff)
            multi[domain][molregno]['target'].append(target)
            multi[domain][molregno]['actId'].append(actId)

  outfile = open('data/multiDict_pKi%s_%s.pkl'%(int(threshold), release),'w')
  pickle.dump(multi, outfile)
  outfile.close()
  return multi          
      
    


   
                                                                         
                                                                         
    
