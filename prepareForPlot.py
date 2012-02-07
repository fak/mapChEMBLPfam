"""
  Function: prepareForPlot(lengths, counts, ran)
  Puts the data into the right format for a boxplot. Needed is a list of lists.
  --------------------
  momo.sander@googlemail.com
"""
def prepareForPlot(lengths, counts, ran):
  import scipy.stats as stat
  data =  {}
  i = 0
  for count in counts:
    try:
      data[count].append(lengths[i])
    except KeyError:
      data[count] = []
      data[count].append(lengths[i])    
    i+=1
    
  ll =[]
  for i in range(ran):
    try:
      ll.append(data[i])
    except KeyError:
      pass
  
  
  for i in range(ran-1):
    x = ll[i]
    y = ll[i+1]
    print "T, P-Value and lower domain count",stat.ttest_ind(x,y, axis = 0), i

  return ll



"""
  Function: getCounts ..gets the domain counts from a specified dictionary 
  pfamDict[target]['count'/'countUnique'] and targets list.  
  --------------------
  momo.sander@googlemail.com
"""                                                  
def getCounts(targets, pfamDict, release, user, pword, host, port):
  import queryDevice
  uniqDomainCount = []
  domainCount= []
  protLength  =[]
  notHuman = []
  lenRatios = []
  for target in targets.keys():
    try:
      pfamDict[target]
    except KeyError:
      print 'not found...', target
      continue
      
    uniqDomainCount.append(pfamDict[target]['countUnique'])
    domainCount.append(pfamDict[target]['count'])
    seq = targets[target]
    print target, seq
    protLength.append(len(seq)-1) 
    ratio = pfamDict[target]['ratio']
    if ratio > 1:
      ratio =1
    lenRatios.append(ratio)
  return (uniqDomainCount, domainCount, protLength, lenRatios)

   
