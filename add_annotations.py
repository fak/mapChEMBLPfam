import os
import queryDevice

####
#### Load parameters.
####

import yaml
# Read config file.
paramFile = open('mpf.yaml')
params = yaml.safe_load(paramFile)
user = params['user']
pword = params['pword']
host = params['host']
port = params['port']
th = params['threshold']
release = params['release']

def add_target_class(level, key, path):
    '''Add target class to each row in a table using a Uniprot accession for release < ChEMBL 15.
    Inputs:
    level: the target class annotation level
    key: the keyword in the header indicating a Uniprot Id
    path: path to the file.
    --------------------
    Felix Kruger
    momo.sander@googlemail.com
    '''
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    out = open('_'.join([path,"sed"]) ,'w')
    out.write('%s\ttarget_class_%s\n'%(lines[0].rstrip('\n'), level))
    header = lines[0].split('\t')
    accessions = {}
    for i, col in enumerate(header):
        if col == key:
            idx = i
            break
    for line in lines[1:]:
        elements = line.split('\t')
        acc = elements[idx]
        accessions[acc] = 0
    accStr = "','".join(map(str, accessions.keys()))
    print len(accessions.keys())
    data = queryDevice.queryDevice("""SELECT protein_accession, %s 
                FROM target_class tc 
                JOIN target_dictionary td 
                  ON td.tid = tc.tid  
                WHERE protein_accession IN('%s')"""%(level, accStr), release, user, pword, host, port)
    for tup in data:
        acc = tup[0]
        targetClass = tup[1]
        accessions[acc] = targetClass
    for line in lines[1:]:
        elements = line.split('\t')
        acc = elements[idx]
        targetClass = accessions[acc]
        out.write("%s\t%s\n"%(line.rstrip('\n'), targetClass ))
    out.close()
    os.system('mv %s %s'% ('_'.join([path,"sed"]), path))


def add_species(key, path):
    '''Add species annotation to each row in a table using a Uniprot accession for release < ChEMBL 15.
    Inputs:
    key: the keyword in the header indicating a Uniprot Id
    path: path to the file.
    --------------------
    Felix Kruger
    momo.sander@googlemail.com
    '''
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    out = open('_'.join([path,"sed"]) ,'w')
    out.write('%s\tspecies_%s\n'%(lines[0].rstrip('\n'), key))
    header = lines[0].split('\t')
    accessions = {}
    for i, col in enumerate(header):
        if col == key:
            idx = i
            break
    for line in lines[1:]:
        elements = line.split('\t')
        acc = elements[idx]
        accessions[acc] = 0
    accStr = "','".join(map(str, accessions.keys()))
    data = queryDevice.queryDevice("""SELECT protein_accession, organism 
		FROM target_dictionary 
		WHERE protein_accession IN('%s') """% accStr, release, user, pword, host, port)
    for tup in data:
        acc = tup[0]
        org = tup[1]
        accessions[acc] = org
    for line in lines[1:]:
        elements = line.split('\t')
        acc = elements[idx]
        org = accessions[acc]
        out.write("%s\t\"%s\"\n"%(line.rstrip('\n'), org))
    out.close()
    os.system('mv %s %s'% ('_'.join([path,"sed"]), path))







