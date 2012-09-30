"""
Function:  master 

Goes through all necessary steps.
  --------------------
  Author:
  Felix Kruger
  fkrueger@ebi.ac.uk
"""  

def master(release, user, pword, host, port): 
  

  ## Get all human protein coding genes from ensembl and a count of all uniqe domains.
  import os
  import pfamDomains
  import mapPfamDomains
  import pdbChembl
  import uniprotChembl
  import analysis

  #Set the threshold for ligand binding to 50 micromolar
  th = 50
  
  # Get Uniprot identifiers for all human proteins.
  os.system("R CMD BATCH --vanilla queryBioMaRt.R")
  # Map Pfam domains and positions to all Uniprot identifiers.
  pfamDomains.pfamDomains(release, user, pword, host, port)
  # Map small molecule binding to PFam domains.
  mapPfamDomains.mapPDs(th, release, user, pword, host, port)
  # Get all ChEMBL interactions in PDB and binding site residues.
  pdbDict = pdbChembl.query(release, user, pword, host, port)
  # Get all ChEMBL interactions in Uniprot and binding site annotation.
  uniprotDict = uniprotChembl.query(release, user, pword, host, port)
  # Analyze the data.
  analysis.analysis(th, release, user, pword, host, port)

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
