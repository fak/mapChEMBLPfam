
"""
  Function:  parse
  --------------------
  assembly of functions to parse text files
  
  momo.sander@googlemail.com
"""                              


def col2keys(path,col, header): 

  i = 0
  if header == True:
    i =1
  infile  = open(path, 'r')
  lines = infile.readlines()
  keyDict = {}
  for line in lines[i:]:
    x = line.split('\t')
    x = x[col].rstrip('\n')  
    keyDict[x] = 0
  return keyDict

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


def rdstatLogs(path): 

  infile = open(path, 'r')
  lines = infile.readlines()
  elements = lines[0].split('\t') 
  al = float(elements[2])
  minx = float(elements[4])
  return(al, minx)



def parse2col(path, header, keyIndex, valIndex):
  dctn = {} 
  i = 0
  if header == True:
    i =1
  infile = open(path, 'r')
  lines = infile.readlines()

  for line in lines[i:]:
    elements = line.split('\t')
    key = elements[keyIndex].rstrip('\n')
    try:
      value = int(elements[valIndex])
    except ValueError:
      value = elements[valIndex].rstrip('\n')
    dctn[key] = value
  return dctn


