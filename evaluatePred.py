"""
  Function:  pdbe
  --------------------
  carries out validation of the algorithm against PDBe.
  
  momo.sander@ebi.ac.uk
"""                                       
def pdbe(pdbDict, keyword,   release): 
                           

  import numpy as np
  out = open('data/%s_PDB_%s.tab'% (keyword, release) , 'w')
  out.write('target\tcmpd\tpPfam\n')
  #humanTargets = humanProtCodUniq.keys()
  
  for target in pdbDict.keys():
    Ts = 0
    Fs = 0
    for cmpdId in pdbDict[target].keys():
      try:
        pdbDict[target][cmpdId][keyword]
      except KeyError:
        #print 'no predictions made for:', target, cmpdId
        continue        
      for pred in pdbDict[target][cmpdId][keyword]:
        if pred:
          Ts +=1
        else:
          Fs +=1
    tot = Ts+Fs
    if tot > 0:
      out.write('%s\t%s\t%s\n'%(target,cmpdId, np.true_divide(Ts, tot))) 
  out.close() 

  return 


def pdbePredicted(pdbDict, keyword,   release): 
                           

  import numpy as np
  out = open('data/%s_PDB_%s.tab'% (keyword, release) , 'w')
  out.write('target\tcmpd\tpPfam\n')
  #humanTargets = humanProtCodUniq.keys()
  
  for target in pdbDict.keys():
    Ts = 0
    Fs = 0
    for cmpdId in pdbDict[target].keys():
      try:
        mapType = pdbDict[target][cmpdId]['mapType']
      except KeyError:
        #print 'no predictions made for:', target, cmpdId
        continue        
      for pred in pdbDict[target][cmpdId][keyword]:
        if pred:
          Ts +=1
        else:
          Fs +=1
    tot = Ts+Fs
    if tot > 0:
      out.write('%s\t%s\t%s\t%s\n'%(target,cmpdId, np.true_divide(Ts, tot), mapType)) 
  out.close() 

  return 
 

"""
  Function:  uniprot
  --------------------
  Carries out validation of the algorithm against Uniprot.
  
  momo.sander@ebi.ac.uk
"""   
def uniprot(bsDict, keyword,  release): 
  import numpy as np
  out = open('data/%s_Uni_%s.tab'% (keyword, release), 'w')
  out.write('target\tpPfam\n')
  for target in bsDict.keys():
    Ts = 0
    Fs = 0
    try:
      bsDict[target][keyword]
    except KeyError:
      #print 'no observations made for:', target
      continue        
    for pred in bsDict[target][keyword]:
      if pred:
        Ts +=1
      else:
        Fs +=1
    tot = Ts+Fs
    if  tot > 0: 
      out.write('%s\t%s\n'%(target, np.true_divide(Ts, tot)))

  out.close()
  return 
  
  
  
"""
  Function:  uniprot
  --------------------
  Carries out validation of the algorithm against Uniprot.
  
  momo.sander@ebi.ac.uk
"""   
def uniprotPredicted(bsDict, keyword,  release): 
  import numpy as np
  out = open('data/%s_Uni_%s.tab'% (keyword, release), 'w')
  out.write('target\tpPfam\n')
  for target in bsDict.keys():
    Ts = 0
    Fs = 0
    try:
      mapType = bsDict[target]['mapType']
    except KeyError:
      #print 'no observations made for:', target
      continue        
    for pred in bsDict[target][keyword]:
      if pred:
        Ts +=1
      else:
        Fs +=1
    tot = Ts+Fs
    if  tot > 0: 
      out.write('%s\t%s\t%s\n'%(target, np.true_divide(Ts, tot), mapType))

  out.close()
  return 
