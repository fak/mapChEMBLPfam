"""
Function:  pfamDomains 

creates pfamDict and pfam_domains
    --------------------
    Author:
    Felix Kruger
    momo.sander@googlemail.com
"""  

def getAllTargets(humanTargets, chemblTargets):
  tDict = {}
 
  for target in chemblTargets:
    tDict[target] = 'chembl'
  for target in humanTargets:
    tDict[target] = 'human'
 
  return tDict



def parse2col(path, header, keyIndex, valIndex):
    dctn = {}
    i = 0
    if header == True:
        i =1
    infile = open(path, 'r')
    lines = infile.readlines()
    for line in lines[i:]:
        elements = line.split('\t')
        key = elements[keyIndex].rstrip('\n')
        try:
            value = int(elements[valIndex])
        except ValueError:
            value = elements[valIndex].rstrip('\n')
        dctn[key] = value
    return dctn


def getUniprotTargets(release, user, pword, host, port):
    import queryDevice
    rawtargets = queryDevice.queryDevice("""SELECT cs.accession, cs.component_id, tid
        FROM component_sequences cs 
            JOIN target_components tc 
            ON tc.component_id = cs.component_id  
        WHERE db_source IN('SWISS-PROT', 'TREMBL')""", release, user, pword, host, port)
    targets= []
    tids = []
    for target in rawtargets:
        targets.append(target[0])
    return targets


def pfamDomains(release, user, pword, host, port): 
  
    import getUniprotTargets
    import getAllTargets
    import getPfamDomains
    import export

    ## Get all ChEMBL targets with a Uniprot accession.
    chemblTargets = getUniprotTargets(release, user, pword, host, port)
    
    ## Read all human protein coding gene names.
    humProtCod = parse2col('data/proteinCoding.tab', True, 1, 0)
    humanTargets = []
    for tstr in humProtCod.keys():
        humanTargets.append(tstr.split(';')[0])
    print "We are dealing with %s human proteins" %len(humanTargets)
    
    ## Generate a list of all targets that are to be fed into the getPfamDomain procedure.
    allTargets = getAllTargets.getAllTargets(humanTargets, chemblTargets)
    allTargets = allTargets.keys()

    ## Get the domains by parsing Pfam. This step takes long and therefore pickles out the domainDict.
    pfamDict = getPfamDomains.getDomains(allTargets, release)  
    
    ## Export the PfamDict as a mysql table.
    export.exportPfamDict(chemblTargets, pfamDict, release, user, pword, host, port)

