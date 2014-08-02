"""
Function:  master 

Goes through all necessary steps.
    --------------------
    Author:
    Felix Kruger
    fkrueger@ebi.ac.uk
"""  

def master(release): 
    

    ## Get all human protein coding genes from ensembl and a count of all uniqe domains.
    import os
    import pfamDomains
    import mapPfamDomains
    import pdbChembl
    import uniprotChembl
    import analysis
    import yaml
    # Read config file.
    paramFile = open('mpf.yaml')
    params = yaml.safe_load(paramFile)
    user = params['user']
    pword = params['pword']
    host = params['host']
    port = params['port']
    th = params['threshold']


    
    # Get Uniprot identifiers for all human proteins.
    #os.system("R CMD BATCH --vanilla queryBioMaRt.R")
    # Map Pfam domains and positions to all Uniprot identifiers.
    #pfamDomains.pfamDomains(release, user, pword, host, port)
    # Map small molecule binding to Pfam domains.
    #mapPfamDomains.mapPDs(th, release, user, pword, host, port)
    # Get all ChEMBL interactions in PDB and binding site residues.
    pdbDict = pdbChembl.query(release, user, pword, host, port)
    # Get all ChEMBL interactions in Uniprot and binding site annotation.
    #uniprotDict = uniprotChembl.query(release, user, pword, host, port)
    # Analyze the data.
    #analysis.analysis(release)

if __name__ == '__main__':
    import sys

    if len(sys.argv) != 2:  # the program name and the two arguments

        sys.exit("Must specify release")

        
    release = sys.argv[1]

    master(release) 
