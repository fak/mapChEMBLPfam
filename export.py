"""
Function:  exportPfamDict(targets,release, pfamDict) 

Exports the PfamDict to the MySQL data base.

momo.sander@googlemail.com
  --------------------

"""

def exportMapsMySQL(propDict, release, user, pword, host, port): 
  import os 
  
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
  out.close()
                                                                               
  os.system("mysql -u%s -p%s -h%s -P%s  -e 'DROP TABLE chembl_%s.map_pfam'"%(user, pword, host, port, release))
  
  os.system("mysql -u%s -p%s -h%s -P%s -e 'CREATE TABLE map_pfam(activity_id VARCHAR(20), domain VARCHAR(100), molregno INT, protein_accession VARCHAR(20), mapType VARCHAR(20))' chembl_%s" %(user, pword, host, port, release))

  os.system("mysqlimport -u%s -p%s -h%s -P%s --lines-terminated-by='\n' --local chembl_%s data/map_pfam.txt"%(user, pword, host, port, release)) 



def exportConflsMySQL(conflicts, release, user, pword, host, port): 
  import os 
  
  out = open('data/conflicts.txt', 'w')
  for conflict in conflicts:
    for target in conflicts[conflict]:
      print target
      out.write('%s\t%s\n'%(target, conflict))
  out.close()
                                                                               
  os.system("mysql -u%s -p%s -h%s -P%s -e 'DROP TABLE chembl_%s.conflicts'"%(user, pword, host, port, release))
  
  os.system("mysql -u%s -p%s -h%s -P%s -e 'CREATE TABLE conflicts (protein_accession "\
            "VARCHAR(20), conflict VARCHAR(200) )' chembl_%s"%(user, pword, host, port, release))

  os.system("mysqlimport -u%s -p%s -h%s -P%s "\
            "--lines-terminated-by='\n' --local chembl_%s data/conflicts.txt"%(user, pword, host, port, release)) 



def exportPfamDict(targets,pfamDict, release, user, pword, host, port):
  import os
  
  print "Writing Pfam Dict to MySQL chembl_%s.pfam_domains"%release
  out = open('data/pfam_domains_%s.txt'%release,'w')
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
  
  os.system("mysql -u%s -p%s -h%s -P%s -e 'DROP TABLE chembl_%s.pfam_domains'" % (user, pword, host, port, release))
  
  os.system("mysql -u%s -p%s -h%s -P%s -e 'CREATE TABLE pfam_domains(protein_accession VARCHAR(20), domain VARCHAR(100),start INT, end INT)' chembl_%s"%(user, pword, host, port, release))

  os.system("cp data/pfam_domains_%s.txt data/pfam_domains.txt" % release)

  os.system("mysqlimport -u%s -p%s -h%s -P%s --lines-terminated-by='\n' --local chembl_%s data/pfam_domains.txt"%(user, pword, host, port ,release))




"""
  Function:  exportProps
  --------------------
  writes the contents of the property dictionary to a table
  
  momo.sander@ebi.ac.uk
"""                                       
def exportProps(selected, threshold, release, user, pword, host, port): 

  import os  
  import queryDevice
  ### Write output to a table.
  out = open('data/cmpdPropssed.tab', 'w')
  out.write('domain\tmolregno\tmolweight\tlogP\tHBA\tHBD\tPSA\tMAPTYPE\n')
  for domain in selected:
    data = queryDevice.queryDevice("SELECT act.molregno, mw_freebase, acd_logp, HBA, HBD, PSA, mpf.mapType FROM compound_properties cp JOIN molecule_dictionary md ON cp.molregno = md.molregno JOIN activities act ON md.molregno = act.molregno JOIN map_pfam mpf ON act.activity_id =mpf.activity_id WHERE mpf.domain ='%s' AND molecule_type = 'Small molecule'"%domain, release, user, pword, host, port) 
    for tup in data:
      out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(domain, tup[0], tup[1], tup[2], tup[3], tup[4], tup[5], tup[6]))
      
  out.close()

  os.system("sed \"s/None/NA/g\" data/cmpdPropssed.tab > \
    data/cmpdProps_pKi%s_chembl%s.tab" %(int(threshold), release))







"""
  Export a table of conflicts for classification by Samuel
"""
def conflicts4Sam():
  import queryDevice
  conflicts = queryDevice.queryDevice("SELECT DISTINCT  mpf.molregno, mpf.domain, dcs.pubmed_id,  con.conflict, mpf.activity_id FROM map_pfam mpf JOIN activities act ON mpf.activity_id = act.activity_id JOIN docs dcs ON dcs.doc_id = act.doc_id JOIN conflicts con ON mpf.protein_accession = con.protein_accession WHERE mapType = 'conflict'", release) 
  confLkp = {}
  for conflict in conflicts:
    confStr = conflict[3]
    confLkp[confStr] = 0 
  for confStr in confLkp.keys():
    out = open('data/forSam_%s.tab'%confStr, 'w')
    out.write('molregno\tpubmed\tprediction\tactivity_id\n')
    for conflict in conflicts:
      if confStr == conflict[3]:
        molregno = conflict[0]
        pubmed = conflict[2]
        domain = conflict[1]
        actId = conflict[4]
        out.write('%s\t%s\t%s\t%s\n'%(molregno, pubmed, domain, actId))
    out.close()
  return

  

