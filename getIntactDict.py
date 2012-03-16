
"""
  Function:  makeIntactDictAllPDBs
  --------------------
  uses a ChEMBL query to make a dictionary of all interaction of PDB ligands 
  with ChEMBL targets.
  
  momo.sander@ebi.ac.uk
"""

def getIntacts(uniDict , release, user, pword, host, port):

  import pickle
  import queryDevice

  chemblIds = uniDict.keys()
  idString = "\',\'".join(chemblIds)

  intactPairs = queryDevice.queryDevice("SELECT md.chembl_id, td.protein_accession \
  FROM activities act \
  JOIN assay2target a2t ON act.assay_id = a2t.assay_id \
  JOIN target_dictionary td ON a2t.tid = td.tid \
  JOIN molecule_dictionary md ON act.molregno = md.molregno \
  JOIN compound_properties cp ON md.molregno = cp.molregno \
  WHERE td.protein_accession IS NOT NULL AND cp.mw_freebase <= 1000 AND md.chembl_id IN('%s')" %idString, release, user, pword, host, port)

  intactDict = {}
  molDict = {}

  for pair in intactPairs:
    chembl_id = pair[0]
    target = pair[1]
    try:
      intactDict[target][chembl_id] = {}
    except KeyError:
      intactDict[target] = {}
      intactDict[target][chembl_id] = {}

  return intactDict
