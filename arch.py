"""
    Script:  arch.py
    Identifies entries in pdbDict for which small molecule binding is mediated
    through more than one domain. The main script is arch.
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
min_res = params['min_res']
min_ratio = params['min_ratio']
####
#### Define functions.
####
def exportPrimers(archs, outfile):
    '''Write out identified architectures in markdown tables.
    Inputs:
    archs -- dictionary of identified architectures.
    '''
    out = open('%s.md' % outfile ,'w')
    out.write('|domain combination|ligand|pdb accession|ratios|# ChEMBL targets|\n')
    out.write('|:-----------|:---------|------------:|:-------|:------------:|\n')
    for arch in sorted(archs.keys()):
        els = arch.split(" $$$ ")
        if len(els) < 2:
            continue
        pString = ", ".join(els)
        for i,ligand in enumerate(archs[arch]['cmpdId']):
            pdbs = {}
            for pdb in archs[arch]['pdb'][i]:
                pdbs[pdb] = 0
            ratios = archs[arch]['ratios'][i]
            rstring = ', '.join(['%.2f' % x for x in ratios])
            pdbString = ", ".join(pdbs.keys()[:min(3,len(pdbs.keys()))])
            ntargs = len(archs[arch]['targets'])
            out.write("|%s|%s|%s|%s|%s|\n"%(pString, ligand,  pdbString, rstring, ntargs))
    out.close()
    outF.close()
    out = open('%s.pkl'%outfile, 'w')
    pickle.dump(archs, out)
    out.close()  


def longPfamDict(pfamDict):
    """Convert pfamDict to long format.
    Inputs:
    pfamDict -- dictionary of Pfam-A domains for Uniprot-ids.
    """
    longPD = {}
    for accession in pfamDict.keys():
        longPD[accession] = {}
        for i, domain in enumerate(pfamDict[accession]['domains']):
            start = pfamDict[accession]['start'][i]
            end = pfamDict[accession]['end'][i]
            for j in range(start, end):
                longPD[accession][j] = domain
    return longPD


def findArch(pdbDict, longPD, min_res, min_ratio):
    """Identify architectures that mediate small molecule binding through
    multiple domain types.
    Inputs:
    pdbDict   -- dictionary, of protein-ligand complexes obtained from 
                 PDBeMotif.
    longPD    -- dictionary, pfamDict in long format.
    min_res   -- int, minimum number of residues per domain.
    min_ratio -- float, minimum ratio of domain contribution to ligand binding
                  over all binding residues.
    """
    archs = {}
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
            arch= ' $$$ '.join(tmpDict.keys())
            ratios = tmpDict.values()
            pdbs = {}
            for pdb in pdbDict[target][cmpdId]['pdb']:
                pdbs[pdb] = 0
            try:
                archs[arch]['ratios'].append(ratios)
                archs[arch]['pdb'].append(pdbs.keys())
                archs[arch]['cmpdId'].append(cmpdId)
            except KeyError:
                archs[arch] = {}
                archs[arch]['pdb'] = [pdbs.keys()]
                archs[arch]['cmpdId'] = [cmpdId]
                archs[arch]['ratios'] = [ratios]
    return archs            



def mapTargets(chemblTargets, pfamDict, archs):
    """
    Function:  mapTargets
    Inputs:
    archs       -- dictionary of identified architectures.    
    chemblTargets --
    pfamDict -- 
    --------------------    
    momo.sander@ebi.ac.uk
    """
    for target in chemblTargets:
        try:
            pfamDict[target]
        except KeyError:
            continue
        for arch in archs.keys():
           arch_doms = arch.split(' $$$ ')
           target_doms = {}
           for dom in arch_doms:
               if dom in pfamDict[target]['domains']:
                   target_doms[dom] = 0
               else:
                   continue   
           if len(arch_doms) == len(target_doms):
              try:
                  archs[arch]['targets'][target]=0
              except KeyError: 
                  archs[arch]['targets']={target:0}
    return archs


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
    archs = findArch(pdbDict, longPD, min_res, min_ratio)
    ## Add targets with given architecture.
    archs = mapTargets(chemblTargets, pfamDict, archs)
    ##  Write architechtures to markdown tables.
    exportPrimers(archs, 'data/interface_%s'%release)  


if __name__ == '__main__':
        import sys

        if len(sys.argv) != 1: # the program name and the two arguments
                print("Parameters are read from mpf.yaml!")
                sys.exit("Parameters are read from mpf.yaml!")


        master() 
