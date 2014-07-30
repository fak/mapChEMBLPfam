"""
    Script:  arch.py
    Identifies entries in pdb_d for which small molecule binding is mediated
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
def export_archs(arch_d, outfile):
    '''Write out identified architectures in markdown tables.
    Inputs:
    arch_d -- dictionary of identified architectures.
    '''
    out = open('%s.md' % outfile ,'w')
    out.write('|domain combination|ligand|pdb accession|ratios|# ChEMBL targets|\n')
    out.write('|:-----------|:---------|------------:|:-------|:------------:|\n')
    for arch in sorted(arch_d.keys()):
        els = arch.split(" $$$ ")
        if len(els) < 2:
            continue
        dom_str = ", ".join(els)
        for i,ligand in enumerate(arch_d[arch]['cmpdId']):
            pdbs = {}
            for pdb in arch_d[arch]['pdb'][i]:
                pdbs[pdb] = 0
            ratios = arch_d[arch]['ratios'][i]
            r_str = ', '.join(['%.2f' % x for x in ratios])
            pdb_str = ", ".join(pdbs.keys()[:min(3,len(pdbs.keys()))])
            n_targs = len(arch_d[arch]['targets'])
            out.write("|%s|%s|%s|%s|%s|\n"%(dom_str, ligand,  pdb_str, r_str, n_targs))
    out.close()
    out = open('%s.pkl'%outfile, 'w')
    pickle.dump(arch_d, out)
    out.close()  


def get_long_pfams(pfam_d):
    """Convert pfam_d to long format.
    Inputs:
    pfam_d -- dictionary of Pfam-A domains for Uniprot-ids.
    """
    long_pfam_d = {}
    for accession in pfam_d.keys():
        long_pfam_d[accession] = {}
        for i, domain in enumerate(pfam_d[accession]['domains']):
            start = pfam_d[accession]['start'][i]
            end = pfam_d[accession]['end'][i]
            for j in range(start, end):
                long_pfam_d[accession][j] = domain
    return long_pfam_d


def get_archs(pdb_d, long_pfam_d, min_res, min_ratio):
    """Identify architectures that mediate small molecule binding through
    multiple domain types.
    Inputs:
    pdb_d   -- dictionary, of protein-ligand complexes obtained from 
                 PDBeMotif.
    long_pfam_d    -- dictionary, pfam_d in long format.
    min_res   -- int, minimum number of residues per domain.
    min_ratio -- float, minimum ratio of domain contribution to ligand binding
                  over all binding residues.
    """
    arch_d = {}
    for target in pdb_d.keys():
        domains = {}
        for domain in long_pfam_d[target].values():
            try:
                domains[domain] += 1
            except KeyError:
                domains[domain] = 1
        for cmpdId in pdb_d[target].keys():
            tmpDict = {}
            for domain in domains:
                dom_res = [x for x in long_pfam_d[target].keys() if long_pfam_d[target][x] == domain and x in pdb_d[target][cmpdId]['position']]
                ndom = len(dom_res)
                nall = len(set(pdb_d[target][cmpdId]['position']))
                ratio = np.true_divide(ndom, nall)
                if ndom >= min_res and ratio >= min_ratio:
                    tmpDict[domain] = ratio
            if len(tmpDict.keys()) <2:
                continue
            arch= ' $$$ '.join(tmpDict.keys())
            ratios = tmpDict.values()
            pdbs = {}
            for pdb in pdb_d[target][cmpdId]['pdb']:
                pdbs[pdb] = 0
            try:
                arch_d[arch]['ratios'].append(ratios)
                arch_d[arch]['pdb'].append(pdbs.keys())
                arch_d[arch]['cmpdId'].append(cmpdId)
            except KeyError:
                arch_d[arch] = {}
                arch_d[arch]['pdb'] = [pdbs.keys()]
                arch_d[arch]['cmpdId'] = [cmpdId]
                arch_d[arch]['ratios'] = [ratios]
    return arch_d            



def add_targets(chembl_targets, pfam_d, arch_d):
    """
    Function:  mapTargets
    Inputs:
    arch_d       -- dictionary of identified architectures.    
    chembl_targets --
    pfam_d -- 
    --------------------    
    momo.sander@ebi.ac.uk
    """
    for target in chembl_targets:
        try:
            pfam_d[target]
        except KeyError:
            continue
        for arch in arch_d.keys():
           arch_doms = arch.split(' $$$ ')
           target_doms = {}
           for dom in arch_doms:
               if dom in pfam_d[target]['domains']:
                   target_doms[dom] = 0
               else:
                   continue   
           if len(arch_doms) == len(target_doms):
              try:
                  arch_d[arch]['targets'][target]=0
              except KeyError: 
                  arch_d[arch]['targets']={target:0}
    return arch_d



def master():
    """
    Function:  master
    --------------------    
    momo.sander@googlemail.com
    """
    ## Load the pdb_d.
    infile = open('data/pdbDict_%s.pkl' %release, 'r')
    pdb_d = pickle.load(infile)
    infile.close()  
    ## Load the pfam_d.
    infile = open('data/protCodPfamDict_%s.pkl' %release, 'r')
    pfam_d = pickle.load(infile)
    infile.close()
    ## Load Uniprot targets.
    chembl_targets = getUniprotTargets.getUniprotTargets(release, user, pword, host, port)
    ## Convert pfam_d to long format.
    long_pfam_d = get_long_pfams(pfam_d)
    ## Identify architectures binding sm through multiple domains.
    arch_d = get_archs(pdb_d, long_pfam_d, min_res, min_ratio)
    ## Add targets with given architecture.
    arch_d = add_targets(chembl_targets, pfam_d, arch_d)
    ##  Write architechtures to markdown tables.
    export_archs(arch_d, 'data/interface_%s'%release)
      

if __name__ == '__main__':
        import sys

        if len(sys.argv) != 1: # the program name and the two arguments
            sys.exit("Parameters are read from mpf.yaml!")


        master() 
