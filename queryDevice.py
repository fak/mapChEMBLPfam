"""
  Function:  queryDevice
  --------------------
  query the ChEMBL
  
  momo.sander@googlemail.com
"""  
def queryDevice(sqlQuery,ChEMBL_version, usr, pword): 

  import MySQLdb
  ChEMBL_version
  
  conn = MySQLdb.connect(host="localhost", user= usr, passwd= pword, \
                         db = "chembl_%s"%(usr, pword,ChEMBL_version))

    	
 	
  query = conn.cursor()
  query.execute(sqlQuery)
  queryResults = query.fetchall()  
  conn.close()

  return queryResults


