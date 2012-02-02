"""
Function:  master 

Goes through all necessary steps.
  --------------------

"""  

def master(release, user, pword): 
  

  ## Get all human protein coding genes from ensembl and a count of all uniqe domains.
  
  os.system("R CMD BATCH --vanilla getPfamStats.R")
  pfamDomains.pfamDomains(release, user, pword)
  mapPfamDomains.mapPDs(release, user, pword)
  pdbDict = pdbChembl.query(release, user, pword)
  uniprotDict = uniprotChembl.query(relese, user, pword)


if __name__ == '__main__':
  import sys
  if len(sys.argv) < 4:  # the program name and the two arguments
    sys.exit("Must specify a dictionary of the form pfamDict[target]['count'] and \
  a lsit of targets")
    
  release = sys.argv[1]
  user = sys.argv[2]
  pword = sys.argv[3]

  master(release, user, pword)
  
  
  
