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
  
  ## Set threshold for all calculations.
  import numpy as np
  threshold = -np.log10(50*10**(-6))


  ## Get all ChEMBL targets with a Uniprot accession.
  import getUniprotTargets
  chemblTargets = getUniprotTargets.getUniprotTargets(release, user, pword, host, port)
  
  ## Read all human protein coding genes
  import parse
  humProtCod = parse.parse2col('data/proteinCoding.tab', True, 1, 0)
  #humanTargets = humanProtCodUniq.keys()
  print "We are dealing with %s human proteins" %len(humProtCod.keys())

  ## Get a list of all human (!) ChEMBL targets
  humChembl = {}
  for target in chemblTargets:
    if target in humProtCod.keys():
      humChembl[target] = 0

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

  ## Load the uniprotDict.
  import pickle
  infile  = open('data/bsDictUniprot_chembl%s.pkl'%release, 'r')
  uniprotDict = pickle.load(infile)
  infile.close()
  print 'number of targets with binding site information', len(uniprotDict.keys())


  ## Load the uniDict.
  import parseUniChem
  uniDict = parseUniChem.parse('data/unichemMappings.txt')

  ## Load the propDict.
  import pickle
  infile = open('data/propDict_%s.pkl'% release, 'r')
  propDict = pickle.load(infile)
  infile.close()

  ####
  #### Generate Plots.
  ####

  ## For each target in PfamDict, calculate the ratio of domain over non-domain regions.
  import getRatioUnstruct
  import writeTable
  import os
  pfamDict = getRatioUnstruct.getRatio(pfamDict, humProtCod, release, user, pword, host, port)
  writeTable.writePfam(pfamDict, humProtCod,humChembl, chemblTargets, release)
  os.system('/ebi/research/software/Linux_x86_64/bin/R-2.11.0 CMD BATCH --vanilla -%s plotPfamStat.R' %release) 


  ## Assess small molecule binding within Pfam domains for PDBe entries.
  import matchData
  import evaluatePred 
  pdbDict = matchData.pdbe(pdbDict,pfamDict, release)
  evaluatePred.pdbe(pdbDict, 'within', release)
  os.system('/ebi/research/software/Linux_x86_64/bin/R-2.11.0 CMD BATCH --vanilla -%s  -%s -%s  ecdf.R' % ('within', "PDB" , release))

  
  ## Assess small molecule binding within Pfam domains for Uniprot entries.  
  import matchData
  import evaluatePred  
  uniprotDict = matchData.uniprot(uniprotDict,pfamDict,  release)
  evaluatePred.uniprot(uniprotDict, 'within', release)
  os.system('/ebi/research/software/Linux_x86_64/bin/R-2.11.0 CMD BATCH --vanilla -%s -%s -%s  ecdf.R' % ('within', "Uni" , release))

  ## Make a barplot of the group sizes for single, multi-one-valid, multi-no-valid,
  ## multi-multi-valid.  
  import groupSize
  import os
  groupsAll = groupSize.groupSize(chemblTargets, pfamDict)
  print "all possible groups (single, none, multi, conflict):",groupsAll
  (single, multi, conflict) = groupSize.groupSizeMap(chemblTargets, release, user , pword, host, port)
  print "all covered targets (single, multi, conflict): ", len(single), len(multi), len(conflict)
  (single, multi, conflict) = groupSize.actSizeMap(chemblTargets, release, user , pword, host, port)
  print "all covered targets (single, multi, conflict): ", len(single), len(multi),len(conflict)


  ## Plot the evaluation of the mappings.
  import queryDevice
  import matchData
  import evaluatePred
  import os

  intacts = queryDevice.queryDevice("SELECT mpf.protein_accession,mpf.domain,mpf.molregno, pfd.start, pfd.end, mpf.maptype, md.chembl_id FROM map_pfam mpf JOIN pfam_domains pfd ON pfd.protein_accession = mpf.protein_accession JOIN molecule_dictionary md ON md.molregno = mpf.molregno WHERE mpf.domain = pfd.domain", release, user, pword, host, port)

  # ...against PDBe  
  pdbDict = matchData.pdbePredicted(pdbDict,  intacts, uniDict)
  evaluatePred.pdbePredicted(pdbDict, 'prediction', release)
  os.system('/ebi/research/software/Linux_x86_64/bin/R-2.11.0 CMD BATCH --vanilla -%s -%s -%s  ecdf.R' % ('prediction', 'PDB' , release))
  # ...against uniprot
  uniprotDict = matchData.uniprotPredicted(uniprotDict,  intacts)
  evaluatePred.uniprotPredicted(uniprotDict, 'prediction', release)
  os.system('/ebi/research/software/Linux_x86_64/bin/R-2.11.0 CMD BATCH --vanilla -%s  -%s -%s  ecdf.R' % ('prediction', "Uni" , release))


  ## Map the overlap
  import overlap
  tholds = [50,10,5,1,0.5,0.1,0.05,0.01,0.005,0.001, 0.0005,0.0001, 0.00005,0.000001]
  overlap.overlap(propDict, tholds, release)  


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
  countFreqs.countLigs(humProtCod.keys(), chemblTargets, release ,user, pword, host, port)
  countFreqs.countDoms(humProtCod.keys(), pfamDict)
  filenames = ['genFreq.tab', 'domLigs.tab', 'targLigs.tab']

  for filename in filenames:
    os.system('/ebi/research/software/Linux_x86_64/bin/R-2.11.0 CMD BATCH --vanilla -%s statPowerLaw.R' %filename)
    al, minx = parse.rdstatLogs('data/powerLawLog%s' % filename)
    freqs = parse.col2intlist('data/%s'%filename, 1, True)
    print len(freqs), minx, al, filename, type(freqs), type(freqs[1])
    plplot.plplot(freqs, minx, al, filename)
    plplotRaw.plplotRaw(freqs, filename) 


  ## Make scatterplot of rank ordered gene counts vs. ligand counts.  
  import prepRank
  import plot

  (genRankL, ligRankL, rectBords) = prepRank.prepRank()
  plot.rankPlot(genRankL, ligRankL, rectBords)   


  ## Plot the ligand properties.
  import export
  import os
  selected = ['7tm_1','Pkinase','Pkinase_Tyr','p450','SNF','Trypsin']
  export.exportProps(selected, threshold, release, user, pword, host, port) 

  filename = 'data/cmpdProps_pKi%s_chembl%s.tab'%(int(threshold), release)
  os.system("/ebi/research/software/Linux_x86_64/bin/R-2.11.0 CMD BATCH --vanilla -%s plotDens.R"%filename)

 
