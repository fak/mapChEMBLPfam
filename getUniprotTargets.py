"""
Function: getUniprotTargets

gets all protein accession with db_source swiss-prot or trembl from the Chembl
target dictionary
momo.sander@googlemail.com
"""
import queryDevice

def getUniprotTargets(release, user, pword, host, port):
    release_number = int(release.split('\_')[1])
    if release_number >= 15:
        rawtargets = queryDevice.queryDevice("""SELECT cs.accession
        FROM component_sequences cs 
            JOIN target_components tc 
            ON tc.component_id = cs.component_id  
        WHERE db_source IN('SWISS-PROT', 'TREMBL')""", release, user, pword, host, port)
    else:
        rawtargets = queryDevice.queryDevice("""SELECT protein_accession
        FROM target_dictionary 
        WHERE db_source IN('SWISS-PROT', 'TREMBL')""", release, user, pword, host, port)
    targets= []
    for target in rawtargets:
        targets.append(target[0])
    return targets 

