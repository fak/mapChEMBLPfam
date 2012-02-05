"""
Function:  master 

Goes through all necessary steps.
  --------------------

"""  

def master(release, user, pword, host, port): 
  

  ## Get all human protein coding genes from ensembl and a count of all uniqe domains.
  import os
  import pfamDomains
  import mapPfamDomains
  import pdbChembl
  import uniprotChembl

  #os.system("R CMD BATCH --vanilla queryBioMaRt.R")
  #pfamDomains.pfamDomains(release, user, pword, host, port)
  mapPfamDomains.mapPDs(release, user, pword, host, port)
  #ipdbDict = pdbChembl.query(release, user, pword, host, port)
  #uniprotDict = uniprotChembl.query(release, user, pword, host, port)


if __name__ == '__main__':
  import sys

  if len(sys.argv) < 5:  # the program name and the two arguments

    sys.exit("Must specify release, user, pword, host, port")

    
  release = sys.argv[1]
  user = sys.argv[2]
  pword = sys.argv[3]
  host = sys.argv[4]
  port = int(sys.argv[5])

  master(release, user, pword, host, port)
  
  
  
