"""
    Script:  mandarch.py
    Identifies domains that only occur in multi-domain proteins. The main
    script is mandarch.
    --------------------
    Felix A Kruger
    momo.sander@ebi.ac.uk
"""                       
####
#### import modules.
####
import operator
import numpy as np
import yaml
import pickle
import getUniprotTargets

####
#### Load parameters.
####
paramFile = open('mpf.yaml')
params = yaml.safe_load(paramFile)
USER = params['user']
PWORD = params['pword']
HOST = params['host']
PORT = params['port']
RELEASE = params['release']
TH = params['threshold']
MIN_RES = params['min_res']
MIN_RATIO = params['min_ratio']

####
#### Define functions.
####
def readfile(path, key_name, val_name):
    """Read two columns from a tab-separated file into a dictionary.
    Inputs:
    path -- filepath
    key_name -- name of the column holding the key
    val_name -- name of the column holding the value
    """
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    lkp = {}
    els = lines[0].rstrip().split('\t')
    for i, el in enumerate(els):
        if el == key_name:
            key_idx = i
        if el == val_name:
            val_idx = i
    for line in lines[1:]:
        elements = line.rstrip().split('\t')
        lkp[elements[key_idx]] = elements[val_idx]
    return  lkp


def get_multi_doms(chembl_targets, pfam_d, valid_doms):
    """Find mandatory multi-domain architectures.
    Inputs:
    path -- filepath
    key_name -- name of the column holding the key
    val_name -- name of the column holding the value
    """
    arch_lkp = {}
    dom_lkp = {}
    for target in chembl_targets:
        try:
            doms = pfam_d[target]['domains']
        except KeyError:
            print "No entry in Pfam for: %s"% target
            continue
        inv_doms = [x for x in doms if x not in valid_doms]
        if len(doms) > 1:#len(inv_doms) == len(doms):
            arch = ', '.join(sorted(doms))
            try:
                arch_lkp[arch] += 1
            except KeyError:
                arch_lkp[arch] = 1
            for dom in set(doms):
                try:
                    dom_lkp[dom] += 1
                except KeyError:
                    dom_lkp[dom] = 1
    return(arch_lkp, dom_lkp)




def export_archs(arch_lkp, valid_doms, path):
    '''Write out multi-domain architectures in markdown tables.
    Inputs:
    arch_lkp -- dictionary of multi-domain architectures.
    '''

    sorted_archs = sorted(arch_lkp.iteritems(), key=operator.itemgetter(1), reverse = True)
    out = open('%s.md' % path ,'w')
    out.write('|architecture|count|mapped|\n')
    out.write('|:-----------|:---------|-----:|\n')
    for arch in sorted_archs:
        mapped = False
        doms = str(arch[0]).split(', ')
        if len([x for x in doms if x in valid_doms]) > 0: 
            mapped = True      
        out.write("|%s|%s|%s|\n"%(arch[0], arch[1], mapped))
   


def export_doms(dom_lkp, valid_doms, path):
    '''Write out identified architectures in markdown tables.
    Inputs:
    dom_lkp -- dictionary of domains occuring in multi-domain architectures.
    '''
    sorted_doms = sorted(dom_lkp.iteritems(), key=operator.itemgetter(1), reverse= True)
    out = open('%s.md' % path ,'w')
    out.write('|domain |count| validated|\n')
    out.write('|:-----------|:-----|-------:|\n')
    for dom in sorted_doms:
        mapped = False
        count = dom[1]
        dom = str(dom[0])
        if dom in valid_doms: 
            mapped = True
        out.write("|%s|%s|%s|\n"%(dom, count, mapped))


def master(version):
    """
    Function:  master
    Run through all steps to identify mandatory muli-domain architectures.
    """
    ## Load the pfamDict.
    infile = open('data/protCodPfamDict_%s.pkl' %RELEASE, 'r')
    pfam_d = pickle.load(infile)
    infile.close()
    # Load the list of validated domains.
    valid_dom_d = readfile('data/valid_pfam_v_%(version)s.tab' % locals(), 'pfam_a', 'pfam_a')
    del valid_dom_d['Pkinase_Tyr']
    valid_doms = valid_dom_d.keys()
    ## Load Uniprot targets.
    chembl_targets = getUniprotTargets.getUniprotTargets(RELEASE, USER, PWORD, HOST, PORT)
    ## Add targets with given architecture.
    (arch_lkp, dom_lkp) = get_multi_doms(chembl_targets, pfam_d, valid_doms)
    ##  Write multi-domain architechtures to markdown tables.
    export_archs(arch_lkp, valid_doms, 'data/multi_dom_archs_%s'% RELEASE)  
    ## Write domains from multi-domain architechtures to markdown tables.
    export_doms(dom_lkp, valid_doms, 'data/multi_dom_doms_%s'% RELEASE)

if __name__ == '__main__':
        import sys

        if len(sys.argv) != 2: # the program name and the two arguments
                sys.exit("""Parameters are read from mpf.yaml but must specify
		 	    version for data/valid_pfam_v_%(version)s.tab""")
        version = sys.argv[1]

        master(version) 
