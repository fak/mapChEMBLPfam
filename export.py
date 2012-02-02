"""
Function:  exportPfamDict(targets,release, pfamDict) 

Exports the PfamDict to the MySQL data base.

momo.sander@googlemail.com
  --------------------

"""

def exportMapsMySQL(propDict, release): 
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
                                                                               
  os.system("mysql -uroot -pgaed5oTh  -e 'DROP TABLE chembl_%s.map_pfam'"% release)
  
  os.system("mysql -u%s -p%s -e 'CREATE TABLE map_pfam(activity_id VARCHAR(20),\
             domain VARCHAR(100), molregno INT, protein_accession VARCHAR(20), \
             mapType VARCHAR(20))' chembl_%s"%(user, pword,release))

  os.system("mysqlimport -u%s -p%s --lines-terminated-by='\n' --local \
             chembl_%s data/map_pfam.txt"%(user, pword,release)) 



def exportConflsMySQL(conflicts, release): 
  import os 
  
  out = open('data/conflicts.txt', 'w')
  for conflict in conflicts:
    for target in conflicts[conflict]:
      print target
      out.write('%s\t%s\n'%(target, conflict))
  out.close()
                                                                               
  os.system("mysql -u%s -p%s  -e 'DROP TABLE chembl_%s.conflicts'"\
            %(user, pword,release))
  
  os.system("mysql -u%s -p%s -e 'CREATE TABLE conflicts (protein_accession \
            VARCHAR(20), conflict VARCHAR(200) )' chembl_%s"%(user, pword,release))

  os.system("mysqlimport -u%s -p%s "\
            "--lines-terminated-by='\n' --local chembl_%s data/conflicts.txt"\
            %(user, pword,release)) 



def exportPfamDict(targets,pfamDict, release, user, pword):
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
  
  os.system("mysql -u%s -p%s -e 'DROP TABLE chembl_%s.pfam_domains'"\
    % (user, pword,release))
  
  os.system("mysql -u%s -p%s -e 'CREATE TABLE pfam_domains(protein_accession \
             VARCHAR(20), domain VARCHAR(100),start INT, end INT)' chembl_%s"\
             %(user, pword,release))
  os.system("cp data/pfam_domains_%s.txt data/pfam_domains.txt" % release)
  os.system("mysqlimport -u%s -p%s --lines-terminated-by='\n' --local \
             chembl_%s data/pfam_domains.txt"%(user, pword,release))
  
 
if __name__ == '__main__':

  export(targets,release, pfamDict) 
