"""
  Function:  prepRank
  --------------------
  Prepare the data for the rank plot. 
  
  momo.sander@ebi.ac.uk
"""     

def prepRank():
  import parse
  from operator import itemgetter

  keyIndex = 0
  valIndex = 1
  header  = True
  genFreq  = parse.parse2col('data/genFreq.tab', header, keyIndex, valIndex) 
  ligFreq  = parse.parse2col('data/domLigs.tab', header, keyIndex, valIndex) 
  ligRankTups  = sorted(ligFreq.items(), key = itemgetter(1), reverse = True)
  genRankTups  = sorted(genFreq.items(), key = itemgetter(1), reverse = True)
  i = 1
  tb = 0
  ligRanks = {}
  for tup in ligRankTups:
    dom = tup[0]
    if tb > tup[1]:
      i += 1
    ligRanks[dom] = i
    tb = tup[1]

  i = 1
  tb = 0
  genRanks = {}
  for tup in genRankTups:
    dom = tup[0]
    if tb > tup[1]:
      i += 1
    genRanks[dom] = i
    tb = tup[1]
    
  genRankL = []
  ligRankL = []
  for dom in ligRanks.keys():
    genRankL.append(genRanks[dom]+1)
    ligRankL.append(ligRanks[dom]+1)



  ## Pepare the rectangular limits
  yr2 = []
  for threshold in [8,40,200, 1000]:

    for tup in ligRankTups:
      if tup[1] < threshold:
        print tup[0], ' has less than %s ligands' %threshold
        yr2.append(ligRanks[tup[0]])
        break
  return genRankL, ligRankL, yr2




