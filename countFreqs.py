"""
  Function:  countFreqs
  --------------------
  Count the targets and compounds per domain using the propDict.
  
  momo.sander@ebi.ac.uk
"""                              
def countLigs(protCod,chemblTargets, release, user, pword, host, port):
  import queryDevice
  ligandPerDom = {}
  ligandPerTarget =  {} 
  for humGene in protCod.keys():
    if not humGene in chemblTargets:
      print 'could not find', humGene
      continue
    qRes = queryDevice.queryDevice("SELECT domain, molregno FROM map_pfam WHERE protein_accession = '%s'"%humGene , release, user, pword, host, port)
    for element in qRes:
      domain = element[0]
      lig = element[1]
      try:
        ligandPerDom[domain][lig]=  0 
      except KeyError:
        ligandPerDom[domain] = {}
        ligandPerDom[domain][lig]=  0
      try:
        ligandPerTarget[humGene][lig] = 0 
      except KeyError:
        ligandPerTarget[humGene] = {}
        ligandPerTarget[humGene][lig] = 0

  domOut = open('data/domLigs.tab','w')
  domOut.write('domain\tfreq\n')
  targetOut = open('data/targLigs.tab','w')
  targetOut.write('domain\tfreq\n')
  for domain in ligandPerDom:
    numLigs = len(ligandPerDom[domain].keys())
    domOut.write('%s\t%s\n'%(domain,numLigs))
  domOut.close()
  for target in ligandPerTarget:
    numLigs = len(ligandPerTarget[target].keys())
    targetOut.write('%s\t%s\n'%(target,numLigs))
  targetOut.close()


def countDoms(protCod, pfamDict):
  domCount = {}
  for humGene in protCod.keys():
    try: 
      pfamDict[humGene]
    except KeyError:
      print 'BiomaRt more up to date than pfamDict for: ', humGene
      continue
    for domain in pfamDict[humGene]['domains']:
      try:
        domCount[domain] += 1
      except KeyError:
        domCount[domain] = 1 
  genOut = open('data/genFreq.tab','w')
  genOut.write('domain\tfreq\n')
  for domain in domCount.keys():
    numGens = domCount[domain]
    genOut.write('%s\t%s\n' %(domain, numGens))
  genOut.close()
                                           
