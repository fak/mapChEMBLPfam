"""
  Function:  makeIntactDictAllPDBs
  --------------------
  uses a ChEMBL query to make a dictionary of all interaction of PDB ligands 
  with ChEMBL targets.
  
  momo.sander@ebi.ac.uk
"""         
    
def mkIntactDictAllPDBs(threshold, release, user, pword):

  import pickle
  import queryDevice

  intactPairs = queryDevice.queryDevice("SELECT cr.molregno, compound_key,protein_accession \
  FROM compound_records cr \
  JOIN activities act ON cr.molregno = act.molregno \
  JOIN assay2target a2t ON act.assay_id = a2t.assay_id \
  JOIN target_dictionary td ON a2t.tid = td.tid \
  WHERE src_id = 6 AND protein_accession IS NOT NULL", release, user, pword)
  
  intactDict = {}
  molDict = {}

  for pair in intactPairs:
    molregno = pair[0]
    cmpdId = pair[1]
    target = pair[2]
   
    molDict[molregno] = cmpdId 
    try:
      intactDict[target][molregno] = {}
    except KeyError:
      intactDict[target] = {}
      intactDict[target][molregno] = {}
          
  return (intactDict, molDict)

