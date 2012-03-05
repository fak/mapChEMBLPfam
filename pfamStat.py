"""
  Function:  pfamStat: creates statistics of the occurence of Pfam domains
  in the ChEMBL database as well as the entire human genome
  --------------------
  momo.sander@googlemail.com
"""                                                  
def pfamStat(chemblTargets, humanTargets, pfamDict, release, user, pword, host, port):
  import prepareForPlot
  import getRatioUnstruct
  reload(prepareForPlot)
  import plot
  reload(plot)

  ## Get a list of all human (!) ChEMBL targets
  humChemblTargets = {}
  for target in chemblTargets:
    if target in humanTargets:
      humChemblTargets[target] = 0

  ## For each target in PfamDict, calculate the ratio of domain over non-domain regions.
  pfamDict = getRatioUnstruct.getRatio(pfamDict,humanTargets, release, user, pword, host, port)
  
  ## Get vectors of domain counts and fractions of structured residues for plotting.
  (humanUniqCounts,humanCounts,humanLengths,lenRatiosHuman) = prepareForPlot.getCounts(humanTargets, pfamDict, release, user, pword, host, port)
    
  (chemblUniqCounts, chemblCounts, chemblLengths,lenRatiosChEMBL)  = prepareForPlot.getCounts(humChemblTargets, pfamDict,  release, user, pword, host, port)
  
  ## Plot Histograms of frequencies. 
  plot.plotHist(humanCounts, chemblCounts, '_%s' %release)
  plot.plotHist(humanUniqCounts, chemblUniqCounts, 'uniq_%s' %release)

  ## Get lenghts of non domain regions and a) ratio of structured over unstructured
  ## residues. b) ratio of the shorter over the longer unstructured regions in 
  ## single domain proteins. For the entire proteome...
  datalol = prepareForPlot.prepareForPlot(lenRatiosHuman, humanCounts,6)
  datalol = datalol[1:]
  plot.ratioPlots(datalol, 'ratiosHuman')
  
  ## Get lenghts of non domain regions and a) ratio of structured over unstructured
  ## residues. b) ratio of the shorter over the longer unstructured regions in 
  ## single domain proteins.For ChEMBL targets...   
  datalol = prepareForPlot.prepareForPlot(lenRatiosChEMBL, chemblCounts,6)
  datalol = datalol[1:]
  plot.ratioPlots(datalol, 'ratiosChembl')  
  











#  humanCountsNoZero = []
#  for count in humanCounts:
#    if not count == 0:
#      humanCountsNoZero.append(count
#  plotHist.plotHist(humanCountsNoZero, chemblCounts, 'totNoZero_%s' %release)
#  
