def getLigandsForTarget(target, release, user, pword, host, port):

  import queryDevice
  ligands = queryDevice.queryDevice("SELECT DISTINCT act.molregno, standard_value,\
         standard_type, standard_units, canonical_smiles, act.relation, act.activity_id \
                                                                                \
                                        FROM activities act \
                                        JOIN assay2target a2t \
                                            ON act.assay_id = a2t.assay_id\
                                        JOIN target_dictionary td \
                                            ON a2t.tid = td.tid \
                                        JOIN assays ass \
                                            ON ass.assay_id = act.assay_id \
                                        JOIN compound_structures cs \
                                            ON cs.molregno=act.molregno \
                                                                        \
                                 WHERE td.protein_accession = '%s' \
                                 AND ass.assay_type='B'  \
                                 AND act.relation ='='    \
                                 AND a2t.multi=0  \
                                 AND a2t.complex=0 \
                                 AND a2t.relationship_type = 'D'\
                                 AND act.standard_type IN('Ki','Kd','IC50', \
                                'EC50','-Log Ki','pKd' , 'pA2', 'pI', 'pKa')" \
                                        %target, release, user, pword, host, port)
  return ligands


"""
  Function:  filterForTarget
  --------------------
"""

def filterForTarget(ligands, threshold):

  import pickle
  import math
  import numpy as np
  #print 'filtering ligands at threshold:', threshold
  ligandList = []
  for data in ligands:
    pAfnty=0
    molregno = data[0]
    standardValue = data[1]

    try:
      standardValue = float(standardValue)
    except:
      continue

    standardValue = float(standardValue)
    if standardValue <=0:
      continue
    standardType = data[2]
    standardUnit = data[3]
    smiles = data[4]
    relation = data[5]
    docId = data[6]
    if standardType in ['Ki','IC50', 'EC50'] and standardUnit == 'nM' and relation == '=':
      pAfnty = (-1) * math.log10(standardValue/float(1000000000))
    elif standardType in ['-Log Ki','pKd', 'pA2', 'pI', 'pKi']:
      pAfnty = standardValue

    else:
      continue

    if pAfnty <= threshold:
      continue

    else:
      ligandList.append([smiles,pAfnty,molregno, docId])
  return ligandList





"""
  Function:  singleDomains
  --------------------
  Identify targets with one binding site containing domain and at least one other domain and train a generic model
  momo.sander@ebi.ac.uk
"""                              

def singleDomains(pfamDict, chemblTargets,threshold, release, user, pword, host, port):
  import pickle
 
  single = {}
  for target in chemblTargets:
    if target in pfamDict:
      if len(pfamDict[target]['domains']) == 1:
        domain = pfamDict[target]['domains'][0]
        ligands = getLigandsForTarget(target, release, user, pword, host, port)
        ligands = filterForTarget(ligands, threshold)
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
              single[domain][molregno]['pAfnty']=[aff]
              single[domain][molregno]['target']=[target]
	      single[domain][molregno]['actId'] = [actId]
              single[domain][molregno]['smiles']=smiles
            except KeyError:
              single[domain]={}
              single[domain][molregno] = {}
              single[domain][molregno]['pAfnty']=[aff]
              single[domain][molregno]['target']=[target]
              single[domain][molregno]['actId'] = [actId]
              single[domain][molregno]['smiles']=smiles
  
  outfile = open('data/singleDict_pKi%s_%s.pkl' %(int(threshold), release),'w')
  pickle.dump(single, outfile)
  outfile.close()  

  return single


                                          
                                                                         
    
