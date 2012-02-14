"""
  Function:  queryDevice
  --------------------
  query the ChEMBL
  
  momo.sander@googlemail.com
"""  
def queryDevice(sqlQuery,ChEMBL_version, usr, pword, hostname, portid ): 

  import MySQLdb
  ChEMBL_version
  
  conn = MySQLdb.connect(host= hostname, user= usr, passwd= pword, db = "chembl_%s"% ChEMBL_version, port = portid )

    	
 	
  query = conn.cursor()
  query.execute(sqlQuery)
  queryResults = query.fetchall()  
  conn.close()

  return queryResults


