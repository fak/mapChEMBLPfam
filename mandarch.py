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
import queryDevice
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
def get_el_targets(RELEASE, USER, PWORD, HOST, PORT):
    data = queryDevice.queryDevice("""
            SELECT DISTINCT td.protein_accession, dc.dc, COUNT(DISTINCT activity_id) 
            FROM target_dictionary td
            JOIN(
                      SELECT pd.protein_accession, COUNT(protein_accession) as dc 
                      FROM pfam_domains pd 
                      GROUP BY pd.protein_accession
                     ) as dc
              ON dc.protein_accession = td.protein_accession
            JOIN assay2target a2t 
              ON td.tid = a2t.tid 
            JOIN assays a 
              ON a.assay_id = a2t.assay_id 
            JOIN activities act 
              ON act.assay_id = a.assay_id
            WHERE act.standard_type IN('Ki','Kd','IC50','EC50', 'AC50') 
            AND a2t.multi=0 
            AND a2t.complex =0 
            AND a2t.relationship_type = 'D' 
            AND act.relation = '=' 
            AND assay_type = 'B' 
            AND standard_value < 50000
            AND dc.dc > 1
            GROUP BY td.protein_accession ORDER BY COUNT(DISTINCT activity_id)""", RELEASE, USER, PWORD, HOST, PORT)
    return data



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


def get_multi_doms(el_targets, pfam_d):
    """Find mandatory multi-domain architectures.
    Inputs:
    el_targets -- list of eligible targets
    pfam_d -- dictionary of Pfam domains
    """
    act_lkp = {}
    arch_lkp = {}
    dom_lkp = {}
    for target in el_targets:
        try:
            doms = pfam_d[target[0]]['domains']
        except KeyError:
            print "No entry in Pfam for: %s"% target
            continue
        #inv_doms = [x for x in doms if x not in valid_doms]
        if len(doms) <= 1:
	    print 'impossible!'#len(inv_doms) == len(doms):
        arch = ', '.join(sorted(doms))
        try:
            arch_lkp[arch] += 1
            act_lkp[arch] += target[2]
        except KeyError:
            arch_lkp[arch] = 1
            act_lkp[arch] = target[2]
        for dom in set(doms):
            try:
                dom_lkp[dom] += 1
            except KeyError:
                dom_lkp[dom] = 1
    return(arch_lkp, dom_lkp, act_lkp)


    
def get_doms(el_targets):
    """Find multi-domain architectures.
    Inputs:
    el_targets -- list of eligible targets
    """
    pfam_lkp = {}
    tids = [x[0] for x in el_targets]
    tidstr = "', '".join(str(t) for t in tids)
    data = queryDevice.queryDevice("""
            SELECT tid, domain_name 
            FROM target_components tc
	    JOIN component_domains cd
	      ON cd.component_id = tc.component_id
            JOIN domains d
	      ON d.domain_id = cd.domain_id
            WHERE tc.tid IN('%s')""" %tidstr, RELEASE, USER, PWORD, HOST, PORT)  
    for ent in data:
        tid = ent[0]
        dom = ent[1]
        try:
            pfam_lkp[tid].append(dom)
        except KeyError:
            pfam_lkp[tid] = [dom]
    return pfam_lkp

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
        doms = str(arch[0]).split(', ')
        mapped = ', '.join([x for x in doms if x in valid_doms]) 
        if len(mapped) == 0:
            mapped = False
        out.write("|%s|%s|%s|\n"%(arch[0], arch[1], mapped))


def export_network(arch_lkp, valid_doms, path):
    '''Write out network file.
    Inputs:
    arch_lkp -- dictionary of multi-domain architectures.
    '''
    lkp = {}
    for arch in arch_lkp.keys():
        doms = arch.split(', ')
        count = arch_lkp[arch]
        if type(doms) is str:
            continue
        for i in range(len(doms)-1):
            for j in range(i+1, len(doms)):
                dom_key = ', '.join(sorted([doms[i],doms[j]]))
                try:
                   lkp[dom_key] += count
                except KeyError:
                   lkp[dom_key] = count
    out = open('%s.tab' % path ,'w')
    out.write('dom_1\tdom_2\tcount\n')
    for link in lkp.keys(): 
        doms = str(link).split(', ')
        out.write("%s\t%s\t%s\n"%(doms[0], doms[1], lkp[link]))
    out.close()


def export_attribs(arch_lkp, valid_doms, path):
    '''Write out network file.
    Inputs:
    arch_lkp -- dictionary of multi-domain architectures.
    '''
    out = open('%s.tab' % path ,'w')
    out.write('dom\tvalid\n')
    lkp = {}
    for arch in arch_lkp.keys():
        doms = arch.split(', ')
        for dom in doms:
            valid = False
            if dom in valid_doms:
                valid = True
            lkp[dom] = valid
    for it in lkp.items():
        out.write("%s\t%s\n"%(it[0], it[1]))
    out.close()


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
    ## Load eligible multi-domain targets.
    el_targets = get_el_targets(RELEASE, USER, PWORD, HOST, PORT)
    ## Add targets with given architecture.
    (arch_lkp, dom_lkp, act_lkp) = get_multi_doms(el_targets, pfam_d)
    ##  Write multi-domain architechtures to markdown tables.
    export_archs(arch_lkp, valid_doms, 'data/multi_dom_archs_%s'% RELEASE)  
    ## Write domains from multi-domain architechtures to markdown tables.
    export_doms(dom_lkp, valid_doms, 'data/multi_dom_doms_%s'% RELEASE)
    ## export network file.
    export_network(arch_lkp, valid_doms, 'data/multi_dom_network_%s'% RELEASE)
    ## export network attribute file.
    export_attribs(arch_lkp, valid_doms, 'data/multi_dom_attributes_%s'% RELEASE)
if __name__ == '__main__':
        import sys

        if len(sys.argv) != 2: # the program name and one argument
                sys.exit("""Parameters are read from mpf.yaml but must specify
		 	    version for data/valid_pfam_v_%(version)s.tab""")
        version = sys.argv[1]

        master(version) 
