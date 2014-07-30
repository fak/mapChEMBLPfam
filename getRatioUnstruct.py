"""
  Function:getRatioUnstruct - creates statistics of the occurence of Pfam domains
  in the ChEMBL database as well as the entire human genome
  --------------------
  Author:
  Felix Kruger
  momo.sander@googlemail.com
"""                                                  
def getRatio(pfamDict,humanTargets, release, user, pword, host, port):
 
  import numpy as np
  import queryDevice
  
  for target in pfamDict.keys():
    pfamDict[target]['ratio']='NA'
    try:
      seq = humanTargets[target]
      seq_len  = len(seq)-1
    except KeyError:
      seq= queryDevice.queryDevice("SELECT protein_sequence FROM target_dictionary WHERE protein_accession = '%s'"%target, release, user, pword, host, port)
      try:
        seq_len = len(seq[0][0])-1
      except IndexError:
        continue

    dom_len = 0
    for i in range(len(pfamDict[target]['domains'])):
      start  = pfamDict[target]['start'][i]         
      end  = pfamDict[target]['end'][i]
      ind_dom = end - start
      dom_len += ind_dom
                                      
    ratio = np.true_divide(dom_len,seq_len)
    pfamDict[target]['ratio'] = ratio
           
  return pfamDict
      
    

      
                                            
#  if len(pfamDict[target]['start']) == 1:
#    start  = pfamDict[target]['start'][i]
#    end  = pfamDict[target]['end'][i]
#    pre = start
#    post = length - end
#    unstructuredRatio = np.true_divide(min([pre,post]), max([pre, post]))
#    pfamDict[target]['weighting'] = unstructuredRatio
#  else:
#    pass
                                    
