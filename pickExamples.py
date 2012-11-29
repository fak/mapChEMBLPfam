
"""
Function: pickExamples

 Select predictions for drugs or compounds with research codes.
  --------------
momo.sander@googlemail.com
"""
 

def researchCode(pfamDict, release, user, pword, host, port):
  import time
  import queryDevice
 
  infile = open('data/map_pfam.txt','r')
  lines = infile.readlines()
  
  out = open('data/resCodeExamples.tab', 'w')
  out.write('research_code\tpfam\tuniprot\tnPfam\n')

  for line in lines[1:]:
    elements = line.split('\t')
    molregno = elements[2]
    pfam = elements[1]
    uniprot = elements[3]
    nPfam = len(pfamDict[uniprot]['domains'])
    time.sleep(0.03)
    resCode = queryDevice.queryDevice("SELECT synonyms FROM molecule_synonyms WHERE syn_type = 'RESEARCH_CODE' AND molregno = %s" %molregno, release, user, pword, host, port)
    try:
      resCode = resCode[0][0]
      out.write('%s\t%s\t%s\t%s\n'%(resCode, pfam, uniprot, nPfam))
    except IndexError:
      pass      
  out.close()

  return


def drugs(pfamDict, release, user, pword, host, port):
  import time
  infile = open('data/map_pfam.txt','r')
  lines = infile.readlines()
  
  out = open('data/drugExamples.tab', 'w')
  out.write('ingredient\tpfam\tuniprot\tnPfam\n')
  for line in lines[1:]:
    elements = line.split('\t')
    molregno = elements[2]
    pfam = elements[1]
    uniprot = elements[3]
    nPfam = len(pfamDict[uniprot]['domains'])
    time.sleep(0.03)
    ingredient = queryDevice.queryDevice("SELECT ingredient FROM formulations WHERE  molregno = %s" %molregno, release, user, pword, host, port)
    try:
      ingredient = ingredient[0][0]
      out.write('%s\t%s\t%s\t%s\n'%(ingredient, pfam, uniprot, nPfam))
    except IndexError:
      pass
      
  out.close()




if __name__ == '__main__':
  import sys

  if len(sys.argv) < 5:  # the program name and the two arguments

    sys.exit("Must specify path to file, release, user, pword, host, port")
 
  release = sys.argv[1]
  user = sys.argv[2]
  pword = sys.argv[3]
  host = sys.argv[4]
  port = int(sys.argv[5])

  ## Load the pfamDict.
  import pickle
  inFile = open('data/protCodPfamDict_%s.pkl' %release, 'r')
  pfamDict = pickle.load(inFile)
  inFile.close()

  researchCode(pfamDict, release, user, pword, host, port)

