"""
  Function:  forSam
  --------------------
  Write all conflicts out for Sam to classify.
  
  momo.sander@ebi.ac.uk
""" 



def toSam(conflicts, threshold, user, pword, host, release, port):
  
  import parse
  import getLigands
  import filterForTarget
  import queryDevice 

  conf = {}
  for confStr in conflicts.keys():
    for target in conflicts[confStr]:      
      ligands = getLigands.getLigandsForTarget(target, release, user, pword, host, port)
      ligands = filterForTarget.filterForTarget(ligands, threshold)

      for ligand in ligands:
        molregno = ligand[2] 
        actId = ligand[3]
        pubmed = queryDevice.queryDevice("SELECT pubmed_id FROM docs JOIN activities act ON act.doc_id = docs.doc_id WHERE activity_id = %s"%actId,release, user, pword, host, port)[0][0]
        pubmed = pubmed
        try:
          conf[confStr][molregno]['actId'].append(actId)
          conf[confStr][molregno]['pubmed'].append(pubmed)
          conf[confStr][molregno]['pubmed'] = []
          conf[confStr][molregno]['pubmed'].append(pubmed)
        except KeyError:
          try:
            conf[confStr][molregno] = {}
            conf[confStr][molregno]['actId'] = []
            conf[confStr][molregno]['actId'].append(actId)
            conf[confStr][molregno]['pubmed'] = []
            conf[confStr][molregno]['pubmed'].append(pubmed)
          except KeyError:
            conf[confStr]={}
            conf[confStr][molregno] = {}
            conf[confStr][molregno]['actId'] = []
            conf[confStr][molregno]['actId'].append(actId)
            conf[confStr][molregno]['pubmed'] = []
            conf[confStr][molregno]['pubmed'].append(pubmed)
 
  confLkp = {}
  for confStr in conf.keys():
    confLkp[confStr] = 0
    
  for confStr in confLkp.keys():
    out = open('data/forSam_%s.pred'%confStr, 'w')
    out.write('molregno\tpubmed\tprediction\tactivity_id\n')
    for conflict in conf[confStr]:
      for molregno in conf[confStr].keys():
        pubmed = conf[confStr][molregno]['pubmed']
        domain = 'None'
        actId = conf[confStr][molregno]['actId']
        for act in actId:
          out.write('%s\t%s\t%s\t%s\n'%(molregno, pubmed[0], domain, act))
    out.close()

