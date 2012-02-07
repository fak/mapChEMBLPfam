"""
Function:  pfamDomains 

creates pfamDict and pfam_domains
  --------------------
  Author:
  Felix Kruger
  momo.sander@googlemail.com
"""  

def pfamDomains(release, user, pword, host, port): 
 
  import getUniprotTargets
  import parse
  import getAllTargets
  import getPfamDomains
  import export

  ## Get all ChEMBL targets with a Uniprot accession.
  chemblTargets = getUniprotTargets.getUniprotTargets(release, user, pword, host, port)
  
  ## Read all human protein coding genes
  humanProtCodUniq = parse2col('data/proteinCoding.tab', True, 1, 0)
  humanTargets = humanProtCodUniq.keys()
  print "We are dealing with %s human proteins" %len(humanTargets)
  
  ## Generate a list of all targets that are to be fed into the getPfamDomain procedure.
  allTargets = getAllTargets.getAllTargets(humanTargets, chemblTargets)
  allTargets = allTargets.keys()

  ## Get the domains by parsing Pfam. This step takes long and therefore pickles out the domainDict.
  pfamDict = getPfamDomains.getDomains(allTargets, release)  
  
  ## Export the PfamDict as a mysql table.
  export.exportPfamDict(chemblTargets, pfamDict, release, user, pword, host, port)

