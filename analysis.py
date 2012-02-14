"""
Function:  analysis 

Carries out the analysis of the data generated in the previous steps.
  --------------------
  Author:
  Felix Kruger
  momo.sander@googlemail.com
"""  


def analysis(release, user, pword, host, port):

  ####
  #### Load data.
  ####

  ## Get all ChEMBL targets with a Uniprot accession.
  import getUniprotTargets
  chemblTargets = getUniprotTargets.getUniprotTargets(release, user, pword, host, port)
  
  ## Read all human protein coding genes
  import parse
  humanProtCodUniq = parse.parse2col('data/proteinCoding.tab', True, 1, 0)
  humanTargets = humanProtCodUniq.keys()
  print "We are dealing with %s human proteins" %len(humanTargets)

  ## Load the pfamDict.
  import pickle
  inFile = open('data/protCodPfamDict_%s.pkl' %release, 'r')
  pfamDict = pickle.load(inFile)
  inFile.close() 

  ## Load the pdbDict.
  import pickle
  infile = open('data/pdbDict_chembl%s.pkl' %release, 'r')
  pdbDict = pickle.load(infile)
  infile.close()

  ## Get binding sites for each target.
  import pickle
  infile  = open('data/bsDictUniprot_chembl%s.pkl'%release, 'r')
  uniprotDict = pickle.load(infile)
  infile.close()
  print 'number of targets with binding site information', len(uniprotDict.keys())

  ## Generate a mapping of PDBe identifiers and molregnos.
  import generateMolDict
  molDict = generateMolDict.generateMolDict(release, user, pword, host, port)


  ####
  #### Generate Plots.
  ####


  ## Plot the histogram of domain numbers per protein and the boxplot of the ratios
  ## of structured over unstructured regions.
  import pfamStat
  pfamStat.pfamStat(chemblTargets, humanProtCodUniq, pfamDict, release, user, pword, host, port)

  ## Assess small molecule binding within Pfam domains for PDBe entries.
  import plot
  import matchData
  import evaluatePred 
  pdbDict = matchData.pdbe(pdbDict,pfamDict, release)
  counts = evaluatePred.pdbe(pdbDict,  release)
  plot.plotEmpCDF(counts, 'pdb_chembl%s' %release)
  
  ## Assess small molecule binding within Pfam domains for Uniprot entries.  
  import plot
  import matchData
  import evaluatePred  
  uniProtDict = matchData.uniprot(uniprotDict,pfamDict,  release)
  counts  = evaluatePred.uniprot(uniprotDict, release)
  plot.plotEmpCDF(counts, 'uniprot_chembl%s' %release) 

  ## Make a barplot of the group sizes for single, multi-one-valid, multi-no-valid,
  ## multi-multi-valid.  
  import groupSize
  import os
  groups = groupSize.groupSize(chemblTargets, pfamDict)
  os.system("R CMD BATCH --vanilla -%s barPlot.R"%(groups))

  ## Plot the evaluation of the mappings.
  import queryDevice
  import matchData
  import evaluatePred
  import os
  mapTypes = ['single', 'multi', 'conflict']
  predLs = {}
  for mapType in mapTypes:
    intacts = queryDevice.queryDevice("SELECT mpf.protein_accession,mpf.domain,mpf.molregno,  pfd.start, pfd.end FROM map_pfam mpf JOIN pfam_domains pfd ON pfd.protein_accession = mpf.protein_accession WHERE mpf.maptype = '%s' AND mpf.domain = pfd.domain"% mapType, release, user, pword, host, port)
    predList = matchData.pdbePredicted(pdbDict,  intacts, molDict, release, mapType)
    predLs[mapType] = predList

  specStr = evaluatePred.prepPlot(predLs, mapTypes)
  os.system("R CMD BATCH --vanilla -%s stackBarPlot.R"%(specStr))


  ## Power Law Distribution of domain occurences
  ##  Prepare the data for the power law plot.
  ##  1. Count the targets and compounds per domain using the propDict
  ##  2. Count a human genes per domain using the Pfam dictionary
  ##  3. Plot the power law distributions for all domains and overlay 25 most 
  ##     frequent domains
  import countFreqs
  import plplot
  import plplotRaw
  import parse 
  countFreqs.countLigs(humanTargets, chemblTargets, release ,user, pword, host, port)
  countFreqs.countDoms(humanTargets, pfamDict)
  filenames = ['genFreq.tab', 'domLigs.tab', 'targLigs.tab']

  for filename in filenames:
    os.system('/ebi/research/software/Linux_x86_64/bin/R-2.11.0 CMD BATCH --vanilla -%s statPowerLaw.R' %filename)
    al, minx = parse.rdstatLogs('data/powerLawLog%s' % filename)
    freqs = parse.col2fltlist('data/%s'%filename, 1, True)
    plplot.plplot(freqs, minx, al, filename)
    plplotRaw.plplotRaw(freqs, filename) 


  ## Make scatterplot of rank ordered gene counts vs. ligand counts.  
  import prepRank
  import plot

  (genRankL, ligRankL, rectBords) = prepRank.prepRank()
  plot.rankPlot(genRankL, ligRankL, rectBords)    
