"""
Function:  getBindingSites(targets) 

Takes a list of ChEMBL targets (input) and determines their binding sites from Uniprot.
Output: bsDict

"""
def getBindingSites(targets,release):
  import urllib
  from xml.dom.minidom import parse  
  import xml.dom
  import xml.parsers.expat
  import pickle
  bsDict ={}
  for target in targets[:4]:
    try:
      print "searching binding sites for %s" % target
      bsDict[target] = {}
      bsDict[target]["ligands"] = []
      bsDict[target]["positions"] = []
      opener = urllib.FancyURLopener({})                                     
      f = opener.open("http://www.uniprot.org/uniprot/%s.xml" % target) 
      dom = parse(f)
      for uniprot in dom.childNodes:
        if uniprot.nodeName == 'uniprot':
          for entry in uniprot.childNodes:
            if entry.nodeName == 'entry':
              for feature in entry.childNodes:
                if feature.nodeName == 'feature':
                  if feature.getAttribute('type') == 'binding site':
                    ligand = feature.getAttribute('description') 
                    for location in feature.childNodes:
                      if location.nodeName == 'location':
                        for position in location.childNodes:
                          if position.nodeName == 'position':
                            position = position.getAttribute('position')
                            print target,'\t',ligand,'\t',position
                            bsDict[target]['ligands'].append(ligand)
                            bsDict[target]['positions'].append(position)
      dom.unlink()
    except xml.parsers.expat.ExpatError:
      continue
    if len(bsDict[target]["positions"]) == 0:
      del bsDict[target]
      
      
  out  = open('data/bsDictUniprot_chembl%s.pkl'%release, 'w')
  pickle.dump(bsDict, out)
  out.close()
      
  return bsDict
  

