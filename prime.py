"""
    Script:  prime.py

    Ientifies entries in pdbDict for which small molecule binding is mediated through more than one domain. The main script is prime.

    --------------------
    Felix A Kruger
    momo.sander@ebi.ac.uk
"""                       
####
#### import modules.
####
import numpy as np
import yaml
import pickle
import getUniprotTargets

####
#### Load parameters.
####
paramFile = open('mpf.yaml')
params = yaml.safe_load(paramFile)
user = params['user']
pword = params['pword']
host = params['host']
port = params['port']
release = params['release']
th = params['threshold']

####
#### Define functions.
####
def exportPrimers(primers, outfile):
    '''Write out identified architectures in markdown tables.
    Inputs:
    primers -- dictionary of identified architectures.
    '''
    outF = open('%s_long.md' % outfile ,'w')
    outF.write('|domain combination|pdb accession|ratios|\n')
    outF.write('|:-----------|------------:|:------------:|\n')
    out = open('%s.md' % outfile ,'w')
    out.write('|domain combination|pdb accession|# ChEMBL targets|\n')
    out.write('|:-----------|------------:|:------------:|\n')
    for primer in sorted(primers.keys()):
        els = primer.split(" $$$ ")
        if len(els) < 2:
            continue
        pString = ", ".join(els)
        pdbs = {}
        for i, pdb in enumerate(primers[primer]['pdb']):
            pdbs[pdb[0]] = 0
            ratios = primers[primer]['ratios'][i]
            rstring = ', '.join(['%.2f' % x for x in ratios])
            outF.write("|%s|%s|%s|\n"%(pString, pdb[0], rstring))
        pdbString = ", ".join(pdbs.keys()[:min(3,len(pdbs.keys()))])
        ntargs = len(primers[primer]['targets'])
        out.write("|%s|%s|%s|\n"%(pString, pdbString, ntargs))
    out.close()
    outF.close()
    out = open('%s.pkl'%outfile, 'w')
    pickle.dump(primers, out)
    out.close()  


def longPfamDict(pfamDict):
    '''Convert pfamDict to long format.
    
    Inputs:
    pfamDict -- dictionary of Pfam-A domains for Uniprot-ids.
    '''
    longPD = {}
    for accession in pfamDict.keys():
        longPD[accession] = {}
        for i, domain in enumerate(pfamDict[accession]['domains']):
            start = pfamDict[accession]['start'][i]
            end = pfamDict[accession]['end'][i]
            for j in range(start, end):
                longPD[accession][j] = domain
    return longPD


def findPrimers(pdbDict, longPD, min_res, min_ratio):
    '''Identify architectures that mediate small molecule binding through multiple domain types.
    
    Inputs:
    pdbDict   -- dictionary, of protein-ligand complexes obtained from PDBeMotif.
    longPD    -- dictionary, pfamDict in long format.
    min_res   -- int, minimum number of residues per domain.
    min_ratio -- float, minimum ratio of domain contribution to ligand binding over all binding residues.
    '''
    primers = {}
    for target in pdbDict.keys():
        domains = {}
        for domain in longPD[target].values():
            try:
                domains[domain] += 1
            except KeyError:
                domains[domain] = 1
        for cmpdId in pdbDict[target].keys():
            tmpDict = {}
            for domain in domains:
                dom_res = [x for x in longPD[target].keys() if longPD[target][x] == domain and x in pdbDict[target][cmpdId]['position']]
                ndom = len(dom_res)
                nall = len(set(pdbDict[target][cmpdId]['position']))
                ratio = np.true_divide(ndom, nall)
                if ndom >= min_res and ratio >= min_ratio:
                    tmpDict[domain] = ratio
            if len(tmpDict.keys()) <2:
                continue
            primer= ' $$$ '.join(tmpDict.keys())
            ratios = tmpDict.values()
            pdbs = {}
            for pdb in pdbDict[target][cmpdId]['pdb']:
                pdbs[pdb] = 0
            try:
                primers[primer]['ratios'].append(ratios)
                primers[primer]['pdb'].append(pdbs.keys())
                primers[primer]['cmpdId'].append(cmpdId)
            except KeyError:
                primers[primer] = {}
                primers[primer]['pdb'] = [pdbs.keys()]
                primers[primer]['cmpdId'] = [cmpdId]
                primers[primer]['ratios'] = [ratios]
    return primers            



def mapTargets(chemblTargets, primers):
    """
    Function:  mapTargets
    Inputs:
    primers       -- dictionary of identified architectures.    
    chemblTargets --
    hierDict      --    
    --------------------    
    momo.sander@ebi.ac.uk
    """
    for target in chemblTargets:
        try:
            pfamDict[target]
        except KeyError:
            continue
        for primer in primers.keys():
           primer_doms = primer.split(' $$$ ')
           target_doms = {}
           for dom in primer_doms:
               if dom in pfamDict[target]['domains']:
                   target_doms[dom] = 0
               else:
                   continue   
           if len(primer_doms) == len(target_doms):
              try:
                  primers[primer]['targets'][target]=0
              except KeyError: 
                  primers[primer]['targets']={target:0}
    return primers


def master():
    """
    Function:  master
    
    --------------------    
    momo.sander@ebi.ac.uk
    """
    ## Load the pdbDict.
    infile = open('data/pdbDict_%s.pkl' %release, 'r')
    pdbDict = pickle.load(infile)
    infile.close()  
    ## Load the pfamDict.
    inFile = open('data/protCodPfamDict_%s.pkl' %release, 'r')
    pfamDict = pickle.load(inFile)
    inFile.close()
    ## Load Uniprot targets.
    chemblTargets = getUniprotTargets.getUniprotTargets(release, user, pword, host, port)
    ## Convert pfamDict to long format.
    longPD = longPfamDict(pfamDict)
    ## Identify architectures binding sm through multiple domains.
    primers = findPrimers(pdbDict, longPD, min_res, min_ratio)
    ## Add targets with given architecture.
    primers = mapTargets(chemblTargets, primers)
    ##  Write architechtures to markdown tables.
    exportPrimers(primers, 'data/interface_%s'%release)  


if __name__ == '__main__':
        import sys

        if len(sys.argv) != 1: # the program name and the two arguments
                print("Parameters are read from mpf.yaml!")
                sys.exit("Parameters are read from mpf.yaml!")


        master() 
