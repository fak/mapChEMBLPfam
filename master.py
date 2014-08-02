"""
Function:  master 
    Master script, calls all scripts in correct order to obtain mapping of
    small molecule binding to Pfam-A (and CATH) domains.

    To run, you need access to a  MySQL instance of ChEMBL > chembl_15. 
    --------------------
    Contact:
    Aurelio Moya-Garcia - aurelio.moya@ucl.ac.uk
    Felix Kruger - fkrueger@ebi.ac.uk
"""  
import yaml
import urllib
from xml.dom.minidom import parse  
import xml.dom
import pickle
import os
import queryDevice
import math
import numpy as np
import MySQLdb

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def mysql_query(query, params, config):
    ''' Extract the chembl_id, molecule and type for each compound in a list of 
    molregnos:

    Inputs:
    params          -- tuple containing query parameters
    default-file    -- path to mysql --defaults-file
    database        -- query database

    '''
    mysql = MySQLdb.connect(user=config['user'], passwd=config['pword'], port=config['port'], host=config['host'] )
    c = mysql.cursor()
    c.execute("use {0}".format(config['release']))
    c.execute(query, params)
    q_res = c.fetchall()
    c.close()
    return q_res

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getUniprotTargets(config):
    query = """SELECT cs.accession, cs.component_id, tid
        FROM component_sequences cs 
            JOIN target_components tc 
            ON tc.component_id = cs.component_id  
        WHERE db_source IN('SWISS-PROT', 'TREMBL')"""
    qres = mysql_query(query, (), config)
    targets= []
    for target in qres:
        targets.append(target[0])
    return targets

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getDomains(targets):
    pfamDict ={}
    ## Loop through all targets and get pfam domains.
    errors = []
    for target in targets:
        #print "getting Pfam domains for %s" % target
        pfamDict[target] = {}
        pfamDict[target]["domains"] = []
        pfamDict[target]["start"] = []
        pfamDict[target]["end"] = []
        opener = urllib.FancyURLopener({})                                     
        f = opener.open("http://pfam.sanger.ac.uk/protein/%s?output=xml" % target) 
        dom = parse(f)
        if not dom.getElementsByTagName('sequence'):
            #print "encountered Error for %s" %target
            errors.append(target)
            del pfamDict[target]
            continue
        for pfam in dom.childNodes:
            if pfam.nodeName == 'pfam':
                for entry in pfam.childNodes:
                    if entry.nodeName == 'entry':
                        for matches in entry.childNodes:
                            if matches.nodeName == 'matches':
                                for match in matches.childNodes:
                                    if match.nodeName == 'match':
                                        if match.getAttribute('type') == 'Pfam-A':
                                            pfamDict[target]['domains'].append(match.getAttribute('id'))
                                            for location in match.childNodes:
                                                if location.nodeName == 'location':
                                                    start = location.getAttribute('start')
                                                    end = location.getAttribute('end')
                                                    pfamDict[target]['start'].append(int(start))
                                                    pfamDict[target]['end'].append(int(end))
        dom.unlink()
        # Add domain count.
        pfamDict[target]['count'] = len(pfamDict[target]['domains'])
        # Calculate and add the uniq count of domains. 
        uniqDomains = {}
        for domain in pfamDict[target]['domains']:
            uniqDomains[domain] = 0   
        pfamDict[target]['countUnique'] = len(uniqDomains)
    ## Pickle the PfamDict
    output = open('data/protCodPfamDict_%s.pkl' % config['release'], 'w')
    pickle.dump(pfamDict, output)
    print "encountered Error for", errors
    return pfamDict   

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def exportPfamDict(targets,pfamDict, config):
    print "Writing Pfam Dict to pfam_domains_{0}".format(config['release'])
    out = open('data/pfam_domains_{0}.txt'.format(config['release']),'w')
    for target in targets:
        i =0
        try:
            for domain in pfamDict[target]['domains']:  
                start  = pfamDict[target]['start'][i]
                end = pfamDict[target]['end'][i]
                i += 1
                out.write("%s\t%s\t%s\t%s\n" % (target, domain, start, end)) 
        except KeyError:
            print 'Couldn\'t find ', target,' in PfamDict.' 
    out.close()
    # Deleted the mysql calls, because mysql export is not strictly needed.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def exportMaps(propDict): 
    out = open('data/map_pfam.txt', 'w')
    for domain in propDict.keys():
        for data in propDict[domain]:
            actIds = data[4]
            if actIds != 'manual':
                for i,actId in enumerate(actIds):
                    molregno = data[0]
                    target = data[2][i]
                    mapType = data[5]
                    out.write('%s\t%s\t%s\t%s\t%s\n'%(actId, domain, molregno,target, mapType))
    # Deleted the mysql calls, because mysql export is not strictly needed.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
# Aurelio - I suggest we a query that is not limited to drugs - this way we 
# capture all domains capable of interactions with small molecules. Since 
# newer releases of chembl have the pchembl value, we can actually use that to
# filter, and drop the whole filterForTarget step \o/

def getLigandsForTarget(target, config):  
    query  = """SELECT DISTINCT ms.canonical_smiles, act.pchembl_value, act.molregno, act.activity_id
    FROM component_sequences cs 
    JOIN target_components tc ON cs.component_id=tc.component_id 
    JOIN target_dictionary td ON tc.tid=td.tid 
    JOIN assays ass ON ass.tid = tc.tid
    JOIN activities act  ON ass.assay_id=act.assay_id
    JOIN molecule_dictionary md ON act.molregno=md.molregno
    JOIN compound_structures ms ON act.molregno=ms.molregno 
    WHERE act.pchembl_value >= {threshold}
    AND potential_duplicate IS NULL
    AND(
        data_validity_comment IS NULL
        OR data_validity_comment = 'manually validated'
        )
    AND assay_type = 'B'
    AND relationship_type = 'D'
    AND target_type = 'SINGLE PROTEIN'
    AND standard_relation = '='
    AND accession = %s""".format(threshold=config['threshold'])
    qres = mysql_query(query, (target,), config)
    ligands = []
    for res in qres:
        ligands.append([res[0],res[1],res[2], res[3]])
    return ligands

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def singleDomains(pfamDict, targets, config):
    single = {}
    for target in targets:
        if target in pfamDict:
            if len(pfamDict[target]['domains']) == 1:
                domain = pfamDict[target]['domains'][0]
                ligands = getLigandsForTarget(target, config)
                for ligand in ligands:
                    smiles = ligand[0]
                    aff = ligand[1]
                    molregno = ligand[2]
                    actId = ligand[3]
                    try:
                        single[domain][molregno]['pAfnty'].append(aff)
                        single[domain][molregno]['target'].append(target)
                        single[domain][molregno]['actId'].append(actId)
                        single[domain][molregno]['smiles']=smiles
                    except KeyError:
                        try:
                            single[domain][molregno] = {}
                            single[domain][molregno]['pAfnty']=[aff]
                            single[domain][molregno]['target']=[target]
                            single[domain][molregno]['actId'] = [actId]
                            single[domain][molregno]['smiles']=smiles
                        except KeyError:
                            single[domain]={}
                            single[domain][molregno] = {}
                            single[domain][molregno]['pAfnty']=[aff]
                            single[domain][molregno]['target']=[target]
                            single[domain][molregno]['actId'] = [actId]
                            single[domain][molregno]['smiles']=smiles
    outfile = open('data/singleDict_pKi%s_%s.pkl' %(int(config['threshold']), config['release']),'w')
    pickle.dump(single, outfile)
    outfile.close()  
    return single

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def multiDomain(pfamDict, targets, valid, config):
    multi = {}
    for target in targets:
        candidates = {}
        if target in pfamDict:
            for domain in pfamDict[target]['domains']:
                domain = str(domain)      
                if domain in valid:
                    candidates[domain]=0
        if len(candidates.keys()) == 1 and len(pfamDict[target]['domains']) >1:
            ligands = getLigandsForTarget(target, config)
            domain = candidates.keys()[0]
            for ligand in ligands:
                smiles = ligand[0]
                aff = ligand[1]
                molregno = ligand[2]
                actId = ligand[3]
                try:
                    multi[domain][molregno]['pAfnty'].append(aff)
                    multi[domain][molregno]['target'].append(target)
                    multi[domain][molregno]['actId'].append(actId)
                    multi[domain][molregno]['smiles']=smiles
                except KeyError:
                    try:
                        multi[domain][molregno] = {}
                        multi[domain][molregno]['pAfnty']=[]
                        multi[domain][molregno]['target']=[]
                        multi[domain][molregno]['actId'] = []
                        multi[domain][molregno]['smiles']=smiles
                        multi[domain][molregno]['pAfnty'].append(aff)
                        multi[domain][molregno]['target'].append(target)   
                        multi[domain][molregno]['actId'].append(actId)
                    except KeyError:
                        multi[domain]={}
                        multi[domain][molregno] = {}
                        multi[domain][molregno]['pAfnty']=[]
                        multi[domain][molregno]['target']=[]
                        multi[domain][molregno]['actId'] = []
                        multi[domain][molregno]['smiles']=smiles
                        multi[domain][molregno]['pAfnty'].append(aff)
                        multi[domain][molregno]['target'].append(target)
                        multi[domain][molregno]['actId'].append(actId)
    outfile = open('data/multiDict_pKi%s_%s.pkl'%(int(config['threshold']), config['release']),'w')
    pickle.dump(multi, outfile)
    outfile.close()
    return multi          
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def col2list(path, col, header): 
    i = 0
    if header == True:
        i =1
    infile  = open(path, 'r')
    lines = infile.readlines()
    ll = []
    for line in lines[i:]:
        elements = line.split('\t')
        element = elements[col].rstrip('\n')
        ll.append(element)
    return ll

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def dictionary(dict_x, propDict, blacklist, maptype): 
    ### feed ligDict into blueprint propDict
    for domain in dict_x.keys(): 
        for molregno in dict_x[domain].keys():
            #print 'feeding %s into propDict. Maptype: %s'%(molregno, maptype)       
            medAfnty = np.median(dict_x[domain][molregno]['pAfnty'])
            smiles = dict_x[domain][molregno]['smiles']
            targets = dict_x[domain][molregno]['target']
            actId = dict_x[domain][molregno]['actId']
            lkp = {}
            for i,target in enumerate(targets):
                if target in blacklist:
                    lkp[target] = 0
                    print 'excluding:', target, 'for domain', domain
            for target in lkp.keys():
                del targets[i]
                del actId[i]
            if len(targets) == 0:
                print 'dropping entry: ', molregno
                continue
            try:
                propDict[domain].append([molregno, smiles, targets,medAfnty,actId, maptype])
            except KeyError:
                propDict[domain] = []
                propDict[domain].append([molregno, smiles, targets, medAfnty,actId, maptype])
    return propDict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def addLigs(propDict, maptype, infile):
    infile = open(infile, 'r')
    lines = infile.readlines()
    infile.close()
    molregno = None
    for line in lines:
        elements = line.split('\t')
        domain = elements[0]
        target = elements[1]
        smiles = elements[2]
        aff = elements[3]
        docId = 'manual'
        print 'whitelisting:', domain       
        try:
            propDict[domain].append([molregno, smiles, target, aff, docId, maptype])
            print 'adding %s to %s' %(smiles, domain)
        except KeyError:
            propDict[domain] = []
            propDict[domain].append([molregno, smiles, target, aff,docId, maptype])
    return propDict  

###############################################################################
###############################################################################
def master(): 
    # Read config file.
    configFile = open('mpf.yaml')
    config = yaml.safe_load(configFile)
    configFile.close()
    ## Get all ChEMBL targets with a Uniprot accession.
    targets = getUniprotTargets(config)
    ## Get the domains by parsing Pfam. This step takes long and pickles out the domainDict.
    pfamDict = getDomains(targets)  
    ## Export the PfamDict as a (mysql) table.
    exportPfamDict(targets, pfamDict, config)
    # Get ligands for all targets with a single domain.
    single = singleDomains(pfamDi, config)
    blacklist = col2list('data/blacklist.tab',1, False)  
    propDict = {}
    propDict = dictionary(single, propDict, blacklist, 'single')
    propDict = addLigs(propDict,'manual', 'data/whitelist.tab') 
    ## Extract a list of validated domains.
    valid = propDict.keys() 
    ## Map small molecule binding to multi-domain targets.
    multi = multiDomain(pfamDict, targets, valid, config)
    ## Insert data for multi domain proteins.
    propDict = dictionary(multi, propDict, blacklist, 'multi')
    ## Export the mapping as a (mysql) table.
    outfile = open('data/propDict_%s.pkl' %config['release'], 'w')
    pickle.dump(propDict, outfile)
    exportMaps(propDict)

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 1:  # the program name and the two arguments
        sys.exit("All parameters specified in config file `mpf.yaml`.")
    master() 
