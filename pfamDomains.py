"""
Function:  pfamDomains 

creates pfamDict and pfam_domains
  --------------------

"""  

def pfamDomains(release, user, pword): 
 
  import getUniProtTargets
  import parse
  import getAllTargets
  import getPfamDomains
  import exportPfamDict

  ## Get all ChEMBL targets with a Uniprot accession.
  chemblTargets = getUniProtTargets.getUniprotTargets(release)
  
  ## Read all human protein coding genes
  humanProtCodUniq = parse.col2keys('data/proteinCoding.tab', 0, True)
  humanTargets = humanProtCodUniq.keys()
  print "We are dealing with %s human proteins" %len(humanTargets)
  
  ## Generate a list of all targets that are to be fed into the getPfamDomain procedure.
  allTargets = getAllTargets.getAllTargets(humanTargets)
  allTargets = allTargets.keys()

  ## Get the domains by parsing Pfam. This step takes long and therefore pickles out the domainDict.
  pfamDict = getPfamDomains.getDomains(allTargets, release)  
  
  ## Export the PfamDict as a mysql table.
  export.exportPfamDict(chemblTargets, pfamDict, release, user, pword)

