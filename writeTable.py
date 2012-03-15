"""
  Function:  writeTable
  --------------------
  momo.sander@googlemail.com
"""                                                  
def writePfam(pfamDict,humProtCod, humChembl, chemblTargets, release):

  out = open('data/pfamTable_%s.tab' % release,'w')
  out.write('target\tnDomains\tpPfam\tsource\n')
  for i,source in enumerate([chemblTargets, humProtCod.keys(), humChembl.keys()]):
    for target in source:
      if target not in pfamDict.keys():
        continue
      if pfamDict[target]['ratio'] == 'NA':
        continue
      nDomains = len(pfamDict[target]['domains'])
      pPfam = pfamDict[target]['ratio']
      srcStr = i
      out.write('%s\t%s\t%s\t%s\n'%(target, nDomains, pPfam, srcStr))
  out.close()


