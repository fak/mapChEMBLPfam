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

  os.system("R CMD BATCH --vanilla getPfamStats.R")
  #pfamDomains.pfamDomains(release, user, pword, host, port)
  mapPfamDomains.mapPDs(release, user, pword, host, port)
  pdbDict = pdbChembl.query(release, user, pword, host, port)
  uniprotDict = uniprotChembl.query(release, user, pword, host, port)


if __name__ == '__main__':
  import sys
<<<<<<< HEAD
  if len(sys.argv) < 5:  # the program name and the two arguments
    sys.exit("Must specify release, user, password, host and port for mySQL connection.")
=======
  if len(sys.argv) < 4:  # the program name and the two arguments
    sys.exit("Must specify a dictionary of the form pfamDict[target]['count'] and \
  a lsit of targets")
>>>>>>> 88ce6159722cc0e2584b54dadf026ee7f066fc8a
    
  release = sys.argv[1]
  user = sys.argv[2]
  pword = sys.argv[3]
  host = sys.argv[4]
  port = int(sys.argv[5])

  master(release, user, pword, host, port)
  
  
  
